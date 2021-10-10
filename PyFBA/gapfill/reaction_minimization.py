import copy
import sys
from random import shuffle

import PyFBA
from PyFBA import log_and_message


def accuracy(precision_recall):
    """
    Calculate the accuracy from the precision and recall calculations. See
    https://en.wikipedia.org/wiki/Precision_and_recall

    :param precision_recall: The dict of tp/tn/fp/fn
    :type precision_recall: dict
    :return: The accuracy of the measurement
    :rtype: float
    """

    return 1.0 * (precision_recall['tp'] + precision_recall['tn']) / (sum(list(precision_recall.values())))


def calculate_precision_recall(growth_media, no_growth_media, modeldata, reactions2run, biomass_eqtn):
    """
    Test growth on our positive and negative media. Return the number of positive/negatives that grew.

    :param no_growth_media: Media on which the model should NOT grow
    :type no_growth_media: list of Media sets
    :param growth_media: Media on which the model should grow
    :type growth_media: list of Media sets
    :param modeldata: The model data
    :type modeldata: PyFBA.model_seed.ModelData
    :param reactions2run:The set of reactions to run
    :type reactions2run: set
    :param biomass_eqtn: The biomass equation
    :type biomass_eqtn: PyFBA.metabolism.reaction.Reaction
    :return: A dict of true positives, true negatives, false positives, false negative
    :rtype: dict of str and int
    """
    results = {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0}
    for media in growth_media:
        status, value, growth = PyFBA.fba.run_fba(modeldata, reactions2run, media, biomass_eqtn)
        PyFBA.fba.remove_uptake_and_secretion_reactions(modeldata.reactions)
        if growth:
            results['tp'] += 1
        else:
            results['fn'] += 1

    for media in no_growth_media:
        status, value, growth = PyFBA.fba.run_fba(modeldata, reactions2run, media, biomass_eqtn)
        PyFBA.fba.remove_uptake_and_secretion_reactions(modeldata.reactions)
        if growth:
            results['fp'] += 1
        else:
            results['tn'] += 1

    return results


def iterate_reactions_to_run(base_reactions, optional_reactions, modeldata, media,
                             biomass_eqn, verbose=False):
    """
    Iterate all the elements in optional_reactions and merge them with base reactions, and then test to see which are
    required for growth

    If you provide an empty set for base_reactions, we test every member of optional_reactions to see if it is required
    for growth.

    This will run as many FBA reactions as there are elements in optional_reactions and only test them one at a time.


    :param modeldata: The model data
    :type modeldata: PyFBA.model_seed.ModelData
    :param base_reactions: a set of reactions that are required for the model but that do not result in growth
    :type base_reactions: set
    :param optional_reactions: a set of reactions that when added to the base_reactions set result in
        growth but for which only a subset may or may not be required.
    :type optional_reactions: set
    :param media: our media object
    :type media: set
    :param biomass_eqn: our biomass equation
    :type biomass_eqn: network.reaction.Reaction
    :param verbose: Print more information
    :type verbose: bool
    :return: The list of reactions that need to be added to base_reactions to get growth
    :rtype: list
    """

    num_elements = len(optional_reactions)
    required_optionals = set()
    i = 1

    for r in range(num_elements):
        removed_reaction = optional_reactions.pop()
        r2r = base_reactions.union(optional_reactions).union(required_optionals)
        log_and_message(f"Single reaction iteration {i} of {num_elements}: Attempting without {removed_reaction}: "
                        f"{modeldata.reactions[removed_reaction].equation}", stderr=verbose)
        status, value, growth = PyFBA.fba.run_fba(modeldata, r2r, media, biomass_eqn)
        if not growth:
            log_and_message("Result: REQUIRED", stderr=verbose)
            required_optionals.add(removed_reaction)
        elif verbose:
            log_and_message("Result: NOT REQUIRED", stderr=verbose)
        i += 1

    return list(required_optionals)


def minimize_additional_reactions(base_reactions, optional_reactions, modeldata, media,
                                  biomass_eqn, verbose=False):
    """
    Given two sets, one of base reactions (base_reactions), and one of optional
    reactions we will attempt to minimize the reactions in the optional
    reactions list by repeatedly bisecting the optional set and testing
    the fba. We return a set of the optional reactions that are required
    for the fba to grow.

    :param base_reactions: a set of reactions that are required for the model but that do not result in growth
    :type base_reactions: set
    :param optional_reactions: a set of reactions that when added to the base_reactions set result in
        growth but for which only a subset may or may not be required.
    :type optional_reactions: set
    :param modeldata: the model seed object that includes compounds and reactions
    :type modeldata: PyFBA.model_seed.ModelData
    :param media: our media object
    :type media: set
    :param biomass_eqn: our biomass equation
    :type biomass_eqn: network.reaction.Reaction
    :param verbose: Print more information
    :type verbose: bool
    :return: The set of reactions that need to be added to base_reactions to get growth
    :rtype: set
    """

    base_reactions = set(base_reactions)
    optional_reactions = set(optional_reactions)
    # test that (a) the base_reactions set does not grow and the base_reactions
    # + optional set does grow
    status, value, growth = PyFBA.fba.run_fba(modeldata, base_reactions, media, biomass_eqn)
    if growth:
        log_and_message("The set of 'base' reactions results in growth so we don't need to bisect the optional set",
                        stderr=True)
        return set()

    status, value, growth = PyFBA.fba.run_fba(modeldata, base_reactions.union(optional_reactions), media,
                                              biomass_eqn)
    if not growth:
        raise Exception("'base' union 'optional' reactions does not generate growth. We can not bisect the set\n")

    # first, lets see if we can limit the reactions based on compounds present and still get growth
    limited_rxn = PyFBA.gapfill.limit_reactions_by_compound(modeldata.reactions, base_reactions, optional_reactions)
    status, value, growth = PyFBA.fba.run_fba(modeldata, base_reactions.union(limited_rxn), media,
                                              biomass_eqn)
    if growth:
        if verbose:
            log_and_message(f"Successfully limited the reactions by compound and reduced from "
                            f"{len(optional_reactions)} to {len(limited_rxn)}", stderr=verbose)
        optional_reactions = limited_rxn

    test = True
    tries = 0
    maxtries = 5
    current_rx_list = list(optional_reactions)
    itera = 0
    log_and_message(f"At the beginning the base list has {len(base_reactions)} and"
                    f" the optional list has {len(current_rx_list)} reactions", stderr=verbose)

    left = []
    right = []
    while test:
        itera += 1
        left, right = PyFBA.gapfill.bisections.bisect(current_rx_list)
        # left, right = percent_split(current_rx_list, percent)
        r2r = base_reactions.union(set(left))
        status, value, lgrowth = PyFBA.fba.run_fba(modeldata, r2r, media, biomass_eqn)
        # running the fba takes all the time, so we only run the right half if the left half doesn't grow
        if lgrowth:
            tries = 0
            current_rx_list = left
            if len(left) == 1:
                test = False
                right = []
            log_and_message(f"Iteration: {itera} Try: {tries} Length: {len(left)} and {len(right)} "
                            f"Growth: {lgrowth} and NOT TESTED", stderr=verbose)
        else:
            r2r = base_reactions.union(set(right))
            status, value, rgrowth = PyFBA.fba.run_fba(modeldata, r2r, media, biomass_eqn)
            log_and_message(f"Iteration: {itera} Try: {tries} Length: {len(left)} and {len(right)} "
                            f"Growth: {lgrowth} and {rgrowth}", stderr=verbose)

            if rgrowth:
                # the right list grows so we can use it
                tries = 0
                if len(right) > 5 * len(left) and len(right) != 1:
                    shuffle(right)
                current_rx_list = right
                if len(right) == 1:
                    test = False
                    left = []
            else:
                # neither grows.
                # If there are less than 20 elements we'll just iterate through them all.
                # Otherwise, we can we split the list unevenly and see if we get growth
                uneven_test = True
                if len(current_rx_list) < 20:
                    left = iterate_reactions_to_run(base_reactions, current_rx_list, modeldata, media, biomass_eqn,
                                                    verbose)
                    right = []
                    test = False
                else:
                    percent = 40
                    left, right = PyFBA.gapfill.bisections.percent_split(current_rx_list, percent)
                    while uneven_test and len(left) > 0 and len(right) > 0:
                        # Not testing left side anymore! The left side is always giving some of its
                        # reactions to the right side. Decreasing the number of reactions will never
                        # result in growth if it didn't grow with the larger number of reactions.
                        # Will leave it commented out for now.

                        # r2r = base_reactions.union(set(left))
                        # status, value, lgrowth = PyFBA.fba.run_fba(compounds, reactions, r2r, media, biomass_eqn)
                        r2r = base_reactions.union(set(right))
                        status, value, rgrowth = PyFBA.fba.run_fba(modeldata, r2r, media, biomass_eqn)
                        log_and_message(f"Iteration: {itera} Try: {tries} Length: {len(left)} and {len(right)} "
                                        f"Growth: {lgrowth} and {rgrowth}", stderr=verbose)
                        # if lgrowth:
                        #    tries = 0
                        #    current_rx_list = left
                        #    uneven_test = False
                        # elif rgrowth:
                        if rgrowth:
                            tries = 0
                            current_rx_list = right
                            uneven_test = False
                        else:
                            percent /= 2.0
                            left, right = PyFBA.gapfill.bisections.percent_split(current_rx_list, percent)
                if uneven_test:
                    # we never got growth, so we can't continue
                    # we take another stab and try again
                    tries += 1
                    current_rx_list = left + right
                    shuffle(current_rx_list)
                if tries > maxtries:
                    test = False

    remaining = set(left + right)
    log_and_message(f"There are {len(remaining)} reactions remaining: {remaining}", stderr=verbose)
    return remaining


def minimize_by_accuracy(base_reactions, optional_reactions, modeldata, growth_media, no_growth_media,
                         biomass_eqn, minimum_tp=0, minimum_accuracy=0.50, verbose=False):
    """
    Given two sets, one of base reactions (base_reactions), and one of optional
    reactions we will attempt to minimize the reactions in the optional
    reactions list by repeatedly bisecting the optional set and testing
    the fba. We choose either the set that maximizes the accuracy of the predictions.

     We return a set of the optional reactions that are required
    for the fba to grow.

    :param minimum_tp: Minimum true positives to consider success. If value < 1 we use that as
            the fraction of growth_media conditions that should be used. (e.g. 0.8 -> 80% of len(growth_media))
    :type minimum_tp: float
    :param minimum_accuracy: The minimum accuracy that we will accept as "growth" e.g. to determine when we have reached
     the end!
    :type minimum_accuracy: float
    :param no_growth_media: A list of media conditions where the model does NOT grow
    :type no_growth_media: list of sets
    :param growth_media: A list of media conditions where the model DOES grow
    :type growth_media: list of sets
    :param base_reactions: a set of reactions that are required for the model but that do not result in growth
    :type base_reactions: set
    :param optional_reactions: a set of reactions that when added to the base_reactions set result in
        growth but for which only a subset may or may not be required.
    :type optional_reactions: set
    :param modeldata: The model data
    :type modeldata: PyFBA.model_seed.ModelData
    :param biomass_eqn: our biomass equation
    :type biomass_eqn: network.reaction.Reaction
    :param verbose: Print more information
    :type verbose: bool
    :return: The set of reactions that need to be added to base_reactions to get growth
    :rtype: set
    """

    if minimum_tp < 1:
        minimum_tp *= len(growth_media)

    base_reactions = set(base_reactions)
    optional_reactions = set(optional_reactions)
    # test that (a) the base_reactions set does not grow and the base_reactions
    # + optional set does grow
    base_precision = calculate_precision_recall(growth_media, no_growth_media, modeldata, base_reactions,
                                                biomass_eqn)
    if base_precision['tp'] > minimum_tp:
        msg = f"The set of 'base' reactions results in {base_precision['tp']} positive reactions. " \
              f"Bigger than {minimum_tp} so no need to bisect"
        log_and_message(msg, stderr=verbose)
        return set()

    base_accuracy = accuracy(base_precision)
    if base_accuracy > minimum_accuracy:
        msg = f"The set of 'base' reactions has an accuracy of {base_accuracy}. " \
              f"Bigger than the threshold of {minimum_accuracy}. No need to bisect"
        log_and_message(msg, stderr=verbose)
        return set()

    beginning_precision = calculate_precision_recall(growth_media, no_growth_media, modeldata,
                                                     base_reactions.union(optional_reactions), biomass_eqn)

    beginning_accuracy = accuracy(beginning_precision)

    if beginning_precision['tp'] < minimum_tp:
        msg = f"If we combine all reactions, we have {beginning_precision['tp']} true positives. We can not get more" \
              f" than this, so we can't exceed {minimum_tp}. No point in continuing"
        log_and_message(msg, stderr=verbose)
        return set()

    msg = f"The beginning accuracy is {beginning_accuracy}. We aim to improve this"
    log_and_message(msg, stderr=verbose)

    # first, lets see if we can limit the reactions based on compounds present and get better accuracy
    limited_rxn = PyFBA.gapfill.limit_reactions_by_compound(modeldata.reactions, base_reactions, optional_reactions)
    new_precision = calculate_precision_recall(growth_media, no_growth_media, modeldata,
                                               base_reactions.union(limited_rxn), biomass_eqn)
    new_accuracy = accuracy(new_precision)
    msg = f"The improved accuracy is {new_accuracy}."
    log_and_message(msg, stderr=verbose)

    if new_precision['tp'] > minimum_tp:
        log_and_message(f"Limited reactions by compound from {len(optional_reactions)} to {len(limited_rxn)}",
                        stderr=verbose)
        optional_reactions = limited_rxn

    test = True
    tries = 0
    maxtries = 5
    current_rx_list = list(optional_reactions)
    itera = 0
    msg = f"At the beginning the base list has {len(base_reactions)} and the " \
          f"optional list has {len(current_rx_list)} reactions"
    log_and_message(msg, stderr=verbose)
    left = []
    right = []
    while test:
        itera += 1
        left, right = PyFBA.gapfill.bisections.bisect(current_rx_list)
        log_and_message(f"Lengths: left {len(left)} right {len(right)}", stderr=verbose)
        # left, right = percent_split(current_rx_list, percent)
        r2r = base_reactions.union(set(left))
        l_precision = calculate_precision_recall(growth_media, no_growth_media, modeldata, r2r, biomass_eqn)
        l_accuracy = accuracy(l_precision)

        r2r = base_reactions.union(set(right))
        r_precision = calculate_precision_recall(growth_media, no_growth_media, modeldata, r2r, biomass_eqn)
        r_accuracy = accuracy(r_precision)

        if l_precision['tp'] > minimum_tp and r_precision['tp'] > minimum_tp:
            log_and_message(f"Both left {l_precision['tp']} and right {r_precision['tp']} are above {minimum_tp}",
                            stderr=verbose)
            # now we want to use the one that increases the accuracy the most
            if l_accuracy > r_accuracy:
                current_rx_list = left
            else:
                current_rx_list = right
            tries = 0
        elif l_precision['tp'] > minimum_tp:
            log_and_message(f"Left precision {l_precision['tp']} is above {minimum_tp}", stderr=verbose)
            current_rx_list = left
            tries = 0
        elif r_precision['tp'] > minimum_tp:
            log_and_message(f"Right precision {r_precision['tp']} is above {minimum_tp}", stderr=verbose)
            current_rx_list = right
            tries = 0
        else:
            log_and_message(f"Neither left {l_precision['tp']} nor right { r_precision['tp']} are above {minimum_tp}",
                            stderr=verbose)

            # neither has increased the tp above minimum_tp, so we have split the list too far
            uneven_test = True
            percent = 40
            left, right = PyFBA.gapfill.bisections.percent_split(current_rx_list, percent)
            while uneven_test and len(left) > 0 and len(right) > 0:
                r2r = base_reactions.union(set(left))
                l_precision = calculate_precision_recall(growth_media, no_growth_media, modeldata,
                                                         r2r, biomass_eqn)
                l_accuracy = accuracy(l_precision)

                r2r = base_reactions.union(set(right))
                r_precision = calculate_precision_recall(growth_media, no_growth_media, modeldata,
                                                         r2r, biomass_eqn)
                r_accuracy = accuracy(r_precision)
                msg = f"Iteration: {itera} Try: {tries} Length: {len(left)} and {len(right)} " \
                      f"Growth: {l_precision['tp']} and {r_precision['tp']} Accuracy {l_accuracy} " \
                      f"and {r_accuracy}"
                log_and_message(msg, stderr=verbose)

                if l_precision['tp'] > minimum_tp:
                    tries = 0
                    current_rx_list = left
                    uneven_test = False
                elif r_precision['tp']:
                    tries = 0
                    current_rx_list = right
                    uneven_test = False
                else:
                    percent /= 2.0
                    left, right = PyFBA.gapfill.bisections.percent_split(current_rx_list, percent)

            if uneven_test:
                # we never got acceptable growth, so we can't continue
                # we take another stab and try again
                tries += 1
                current_rx_list = left + right
                shuffle(current_rx_list)
        if tries > maxtries:
            test = False
        if len(current_rx_list) == 1:
            test = False
            left = current_rx_list
            right = []
    remaining = set(left + right)
    log_and_message(f"There are {len(remaining)} reactions remaining: {remaining}", stderr=verbose)
    return remaining


def minimize_reactions(original_reactions_to_run, added_reactions, modeldata, media, biomass_equation, verbose=False):
    """
    Sort thorugh all the added reactions and return a dict of new reactions
    :param original_reactions_to_run: the original set from our genome
    :type original_reactions_to_run: set(PyFBA.metabolism.Reaction)
    :param added_reactions: new reactions we need
    :type added_reactions: list[(str, set(str)]
    :param modeldata: our modeldata object
    :type modeldata: PyFBA.model_seed.ModelData
    :param media: our media object
    :type media: set[PyFBA.metabolism.Compound]
    :param biomass_equation: our biomass equation
    :type biomass_equation: PyFBA.metabolism.Reaction
    :param verbose: more output
    :type verbose: bool
    :return: A dict of the minimal set of reactions and their source
    :rtype: dict[str, str]
    """
    reqd_additional = set()
    print(f"Before we began, we had {len(original_reactions_to_run)} reactions")

    rxn_source = {}

    while added_reactions:
        ori = copy.deepcopy(original_reactions_to_run)
        ori.update(reqd_additional)
        # Test next set of gap-filled reactions
        # Each set is based on a method described above
        how, new = added_reactions.pop()
        sys.stderr.write(f"Testing reactions from {how}\n")

        # Get all the other gap-filled reactions we need to add
        for tple in added_reactions:
            ori.update(tple[1])

        for r in new:
            # remember the source. It doesn't matter if we overwrite, as it will replace later with earlier
            rxn_source[r] = how

        # Use minimization function to determine the minimal
        # set of gap-filled reactions from the current method
        new_essential = minimize_additional_reactions(ori, new, modeldata, media, biomass_equation,
                                                                    verbose=True)
        log_and_message(f"Saved {len(new_essential)} reactions from {how}", stderr=verbose)
        # Record the method used to determine
        # how the reaction was gap-filled
        for new_r in new_essential:
            modeldata.reactions[new_r].is_gapfilled = True
            modeldata.reactions[new_r].gapfill_method = how
        reqd_additional.update(new_essential)

    # add the original set too
    for r in original_reactions_to_run:
        rxn_source[r] = 'genome prediction'

    # Combine old and new reactions and add the source, to return a dict
    return {r: rxn_source[r] for r in original_reactions_to_run.union(reqd_additional)}

