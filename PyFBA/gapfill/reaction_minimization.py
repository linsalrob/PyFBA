import sys
from random import shuffle

import PyFBA


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
        if verbose:
            sys.stderr.write(f"Single reaction iteration {i} of {num_elements}: Attempting without {removed_reaction}: "
                             f"{modeldata.reactions[removed_reaction].equation}\n")
        status, value, growth = PyFBA.fba.run_fba(modeldata, r2r, media, biomass_eqn)
        if not growth:
            if verbose:
                sys.stderr.write("Result: REQUIRED\n")
            required_optionals.add(removed_reaction)
        elif verbose:
            sys.stderr.write("Result: NOT REQUIRED\n")
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
        sys.stderr.write("The set of 'base' reactions results in growth so we don't need to bisect the optional set\n")
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
            sys.stderr.write("Successfully limited the reactions by compound and reduced " +
                             " from {} to {}\n".format(len(optional_reactions), len(limited_rxn)))
        optional_reactions = limited_rxn

    test = True
    tries = 0
    maxtries = 5
    current_rx_list = list(optional_reactions)
    itera = 0
    sys.stderr.write("At the beginning the base list has {} ".format(len(base_reactions)) +
                     " and the optional list has {} reactions\n".format(len(current_rx_list)))

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
            if verbose:
                sys.stderr.write("Iteration: {} Try: {} Length: {} and {}".format(itera, tries, len(left), len(right)) +
                                 " Growth: {} and NOT TESTED\n".format(lgrowth))
        else:
            r2r = base_reactions.union(set(right))
            status, value, rgrowth = PyFBA.fba.run_fba(modeldata, r2r, media, biomass_eqn)
            if verbose:
                sys.stderr.write("Iteration: {} Try: {} Length: {} and {}".format(itera, tries, len(left), len(right)) +
                                 " Growth: {} and {}\n".format(lgrowth, rgrowth))

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
                        if verbose:
                            sys.stderr.write(
                                "Iteration: {} Try: {} Length: {} and {}".format(itera, tries, len(left), len(right)) +
                                " Growth: {} and {}\n".format(lgrowth, rgrowth))
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
    if verbose:
        sys.stderr.write("There are {} reactions remaining: {}\n".format(len(remaining), remaining))
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
        sys.stderr.write("The set of 'base' reactions results in {} ".format(base_precision['tp']))
        sys.stderr.write("positive reactions. Bigger than {} so no need to bisect\n".format(minimum_tp))
        return set()

    base_accuracy = accuracy(base_precision)
    if base_accuracy > minimum_accuracy:
        sys.stderr.write("The set of 'base' reactions has an accuracy of {} ".format(base_accuracy))
        sys.stderr.write("which is bigger than the threshold of {}. No need to bisect\n".format(minimum_accuracy))
        return set()

    beginning_precision = calculate_precision_recall(growth_media, no_growth_media, modeldata,
                                                     base_reactions.union(optional_reactions), biomass_eqn)

    beginning_accuracy = accuracy(beginning_precision)

    if beginning_precision['tp'] < minimum_tp:
        sys.stderr.write(f"If we combine all reactions, we have get {beginning_precision['tp']} true positives\n")
        sys.stderr.write(f"We can not get more than this, so we can't exceed {minimum_tp}\n")
        sys.stderr.write("No point in continuing\n")
        return set()

    sys.stderr.write(f"The beginning accuracy is {beginning_accuracy}. We aim to improve this\n")

    # first, lets see if we can limit the reactions based on compounds present and get better accuracy
    limited_rxn = PyFBA.gapfill.limit_reactions_by_compound(modeldata.reactions, base_reactions, optional_reactions)
    new_precision = calculate_precision_recall(growth_media, no_growth_media, modeldata,
                                               base_reactions.union(limited_rxn), biomass_eqn)
    new_accuracy = accuracy(new_precision)
    sys.stderr.write(f"The improved accuracy is {new_accuracy}.\n")

    if new_precision['tp'] > minimum_tp:
        if verbose:
            sys.stderr.write(f"Limited reactions by compound from {len(optional_reactions)} to {len(limited_rxn)}\n")
        optional_reactions = limited_rxn

    test = True
    tries = 0
    maxtries = 5
    current_rx_list = list(optional_reactions)
    itera = 0
    if verbose:
        sys.stderr.write(f"At the beginning the base list has {len(base_reactions)} " +
                         f" and the optional list has {len(current_rx_list)} reactions\n")
    left = []
    right = []
    while test:
        itera += 1
        left, right = PyFBA.gapfill.bisections.bisect(current_rx_list)
        if verbose:
            sys.stderr.write("Lengths: left {} right {}\n".format(len(left), len(right)))
        # left, right = percent_split(current_rx_list, percent)
        r2r = base_reactions.union(set(left))
        l_precision = calculate_precision_recall(growth_media, no_growth_media, modeldata, r2r, biomass_eqn)
        l_accuracy = accuracy(l_precision)

        r2r = base_reactions.union(set(right))
        r_precision = calculate_precision_recall(growth_media, no_growth_media, modeldata, r2r, biomass_eqn)
        r_accuracy = accuracy(r_precision)

        if l_precision['tp'] > minimum_tp and r_precision['tp'] > minimum_tp:
            if verbose:
                sys.stderr.write("Both left {} and right {} are above {}\n".format(l_precision['tp'],
                                                                                   r_precision['tp'], minimum_tp))
            # now we want to use the one that increases the accuracy the most
            if l_accuracy > r_accuracy:
                current_rx_list = left
            else:
                current_rx_list = right
            tries = 0
        elif l_precision['tp'] > minimum_tp:
            if verbose:
                sys.stderr.write("Left {} is above {}\n".format(l_precision['tp'], minimum_tp))
            current_rx_list = left
            tries = 0
        elif r_precision['tp'] > minimum_tp:
            if verbose:
                sys.stderr.write("Right {} is above {}\n".format(r_precision['tp'], minimum_tp))
            current_rx_list = right
            tries = 0
        else:
            if verbose:
                sys.stderr.write("Neither left {} nor right {} are above {}\n".format(l_precision['tp'],
                                                                                      r_precision['tp'], minimum_tp))
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
                if verbose:
                    sys.stderr.write(
                        f"Iteration: {itera} Try: {tries} Length: {len(left)} and {len(right)}" +
                        f" Growth: {l_precision['tp']} and {r_precision['tp']} " +
                        f"Accuracy {l_accuracy} and {r_accuracy}\n")
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
    if verbose:
        sys.stderr.write("There are {} reactions remaining: {}\n".format(len(remaining), remaining))
    return remaining
