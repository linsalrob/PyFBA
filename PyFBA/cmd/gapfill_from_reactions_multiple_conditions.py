import argparse
import copy
import os
import sys

import PyFBA
from PyFBA import log_and_message

"""
This code is designed to exemplify some of the gap-filling approaches. If you start with an ungapfilled set of
reactions, we iteratively try to build on the model until it is complete, and then we use the bisection code
to trim out reactions that are not necessary, to end up with the smallest set of reactions.

We will handle multiple positive files and negative files and try and resolve in the most meaningful way (hopefully!)

"""

modeldata = PyFBA.model_seed.ModelData()


def minimize_reactions(original_reactions_to_run, added_reactions, modeldata, media, biomass_equation, verbose=False):
    reqd_additional = set()
    print(f"Before we began, we had {len(original_reactions_to_run)} reactions")

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

        # Use minimization function to determine the minimal
        # set of gap-filled reactions from the current method
        new_essential = PyFBA.gapfill.minimize_additional_reactions(ori, new, modeldata, media, biomass_equation,
                                                                    verbose=True)
        log_and_message(f"Saved {len(new_essential)} reactions from {how}", stderr=verbose)
        # Record the method used to determine
        # how the reaction was gap-filled
        for new_r in new_essential:
            modeldata.reactions[new_r].is_gapfilled = True
            modeldata.reactions[new_r].gapfill_method = how
        reqd_additional.update(new_essential)

    # Combine old and new reactions
    return original_reactions_to_run.union(reqd_additional)


def update_r2r(old, new, why, verbose=False):
    """
    Update the reactions to run and log the changes
    :param old: the initial reactions to run
    :param new: the new reactions to add
    :param why: the step we are at
    :param verbose: more output
    :return: a set of reactions to run
    :rtype: set[str]
    """
    before = len(old)
    old.update(new)
    msg = f"Before updating reactions from {why}: {before} reactions, after {len(old)} reactions"
    log_and_message(msg, stderr=verbose)
    return old


def resolve_additional_reactions(ori_reactions, adnl_reactions, md, gm, ngm, biomass_eqn,
                                 minimum_tp=0, minimum_accuracy=0.5, verbose=False):
    """
    Iteratively resolve additional reactions that are required
    :param ori_reactions: the set of original reactions that form the base of the model
    :type ori_reactions: Set[str]
    :param adnl_reactions: a list of tuples of how the reactions were suggested, and the set of additional reactions
    :param md: The modeldata object
    :type md: PyFBA.model_seed.ModelData
    :param gm: A list of media where we should grow
    :type gm: List[Set[PyFBA.metabolism.CompoundWithLocation]]
    :param ngm: A list of media where we should NOT grow
    :type ngm: List[Set[PyFBA.metabolism.CompoundWithLocation]]
    :param biomass_eqn: our biomass object
    :type biomass_eqn: PyFBA.metabolism.Reaction
    :param minimum_tp: minimum true positives allowed
    :type minimum_tp: float
    :param minimum_accuracy:  minimum accuracy desired
    :type minimum_accuracy: float
    :param verbose:  more output
    :return: set of additional reactions from all of the added_reactions
    :rtype: set
    """

    reqd_additional = set()

    while adnl_reactions:
        ori = copy.copy(ori_reactions)
        ori.update(reqd_additional)
        (how, new) = adnl_reactions.pop()
        sys.stderr.write("Testing suggestions from " + how + "\n")
        # get all the other reactions we need to add
        for tple in adnl_reactions:
            ori.update(tple[1])
        new_essential = PyFBA.gapfill.minimize_by_accuracy(ori, new, md, gm, ngm,
                                                           biomass_eqn, minimum_tp, minimum_accuracy, verbose=verbose)
        for new_r in new_essential:
            md.reactions[new_r].is_gapfilled = True
            md.reactions[new_r].gapfill_method = how
        reqd_additional.update(new_essential)

    return reqd_additional


def read_reactions(reaction_file, verbose=False):
    """
    Read the reactions file and return a set of reactions
    :param reaction_file: the reactions file
    :param verbose: more output
    :return: a set of reactions
    :rtype: set[str]
    """

    reactions2run = set()
    with open(reaction_file, 'r') as f:
        for li in f:
            if li.startswith('#'):
                continue
            if "biomass" in li.lower():
                if verbose:
                    sys.stderr.write("Biomass reaction was skipped from the list as it is auto-imported\n")
                continue
            r = li.strip()
            if r in modeldata.reactions:
                reactions2run.add(r)
    return reactions2run


def measure_accuracy(why, growth_media, no_growth_media, reactions, added_reactions,
                     biomass_eqtn, min_growth_conditions, reaction_source, output, verbose=False):
    """
    Measure the accuracy of this set of reactions
    :param output: the name of the output file
    :type output: str
    :param reaction_source: the dict of where the reactions came from
    :type reaction_source: dict[str, str]
    :param min_growth_conditions: The fraction of reactions we should pass
    :type min_growth_conditions: float
    :param why: why are we doing this
    :type why: str
    :param growth_media: list of media where the organism should grow
    :type growth_media: list[set[PyFBA.metabolism.Compound]]
    :param no_growth_media: list of media where the organism should NOT grow
    :type no_growth_media: list[set[PyFBA.metabolism.Compound]]
    :param reactions: the original reactions to run
    :type reactions: set[str]
    :param added_reactions: the added reactions to run
    :type added_reactions: list[(str, set[str])]
    :param biomass_eqtn: the biomass equation
    :type biomass_eqtn: PyFBA.metabolism.Reaction
    :param verbose: more output
    :type verbose: bool
    :return: True if we have succeeded in positive rate exceeeding min_growth_conditions (ie. we should exit)

    """

    global modeldata
    pr = PyFBA.gapfill.calculate_precision_recall(growth_media, no_growth_media, modeldata, reactions, biomass_eqtn)
    acc = PyFBA.gapfill.reaction_minimization.accuracy(pr)

    msg = f"Testing accuracy: {why} "
    for t in ['tp', 'fp', 'tn', 'fn']:
        msg += f"{t}: {pr[t]} "
    msg += f"Accuracy: {acc}\n"
    log_and_message(msg, stderr=verbose)

    if pr['tp'] > min_growth_conditions:
        log_and_message(f"At step {why} {pr['tp']} is bigger than {min_growth_conditions} : we are done!",
                        stderr=verbose)
        # rediscover the original set of reactions. We could also pass this set :)
        all_reactions = {r for ele in added_reactions for r in ele[1]}
        original_reactions = reactions.difference(all_reactions)
        additions = resolve_additional_reactions(original_reactions, added_reactions, modeldata,
                                                 growth_media, no_growth_media, biomass_eqtn,
                                                 minimum_tp=min_growth_conditions, verbose=verbose)

        with open(output, 'w') as out:
            for r in original_reactions.union(additions):
                if r not in reaction_source:
                    reaction_source[r] = "UNKNOWN??"
                out.write(f"{r}\t{reaction_source[r]}\n")
        sys.exit(0)
    return pr['tp']



def multiple_gapfill(reactions, positive, negative, min_growth_conditions, close, genome_type, output, verbose=False):
    """
    Run multiple gap filling operations and try to resove positive/negative growth
    :param output: Output file name
    :type output: str
    :param reactions: A list of reactions we think are in the genome
    :type reactions: set[str]
    :param positive: A list of media names where the genome grows
    :type positive: list[str]
    :param negative: A list of media names where the genome does not grow
    :type negative: list[str]
    :param min_growth_conditions: The minimum number of conditions where we want growth to consider success
    :type min_growth_conditions: float
    :param close: A list of close genomes used for gapfilling
    :type close: list[str]
    :param genome_type: The genome type (gram +ve, -ve, etc)
    :type genome_type: str
    :param verbose: more output
    :type verbose: bool
    :return: The set of reactions that meets the conditions
    :rtype: set[str]
    """

    if close is None:
        close = list()
    global modeldata
    growth_media = [PyFBA.parse.read_media.find_media_file(m, modeldata=modeldata, verbose=verbose) for m in positive]
    no_growth_media = [PyFBA.parse.read_media.find_media_file(m, modeldata=modeldata, verbose=verbose) for m in negative]
    biomass_eqtn = PyFBA.metabolism.biomass.biomass_equation(genome_type)

    min_growth_conditions = 1.0 * min_growth_conditions * len(growth_media)

    log_and_message(f"Looking for growth on at least {min_growth_conditions} media", stderr=verbose)

    reaction_source = {r: 'original_reaction' for r in reactions}

    max_tp = measure_accuracy('Initial run',  growth_media, no_growth_media, reactions, [], biomass_eqtn,
                              min_growth_conditions, reaction_source, output, verbose)
    best_reactions = copy.deepcopy(reactions)

    added_reactions = []

    #############################################################################################
    #                                       Gapfilling                                          #
    #                                                                                           #
    #  We do this in the order:                                                                 #
    #     essential reactions: because you need to have these, but it is stronger evidence if   #
    #              your friends have it too!                                                    #
    #     linked reactions: connections to other reactions in our set                           #
    #     media: because you should be importing everything in the media                        #
    #     closely related organisms: because you should have roles your friends have            #
    #     subsystems: to complete things you already have                                       #
    #     orphans: to make sure everything is produced/consumed                                 #
    #                                                                                           #
    #     At the moment we don't run these approaches, but easy to add back again!              #
    #     probability: because there are other reactions we can add                             #
    #     reactions with proteins: to make sure you can at least grow on the media              #
    #                                                                                           #
    #############################################################################################

    #############################################################################################
    #                                       ESSENTIAL PROTEINS                                  #
    #############################################################################################

    log_and_message("Gap filling from Essential Reactions", stderr=verbose)
    essential_reactions = PyFBA.gapfill.suggest_essential_reactions()
    for r in essential_reactions:
        modeldata.reactions[r].reset_bounds()
    added_reactions.append(("essential", essential_reactions))
    reactions = update_r2r(reactions, essential_reactions, "ESSENTIAL REACTIONS")

    tp = measure_accuracy('Essential reactions',  growth_media, no_growth_media, reactions, added_reactions,
                          biomass_eqtn, min_growth_conditions, reaction_source, output, verbose)
    if tp > max_tp:
        best_reactions = copy.deepcopy(reactions)
        max_tp = tp

    #############################################################################################
    #                                       LINKED REACTIONS                                    #
    #############################################################################################

    log_and_message("Gap filling from Linked Reactions", stderr=verbose)
    linked_reactions = PyFBA.gapfill.suggest_linked_reactions(modeldata, reactions)
    for r in linked_reactions:
        modeldata.reactions[r].reset_bounds()
    added_reactions.append(("linked_reactions", linked_reactions))
    reactions = update_r2r(reactions, linked_reactions, "LINKED REACTIONS")

    tp = measure_accuracy('Linked reactions',  growth_media, no_growth_media, reactions, added_reactions, biomass_eqtn,
                          min_growth_conditions, reaction_source, output, verbose)
    if tp > max_tp:
        best_reactions = copy.deepcopy(reactions)
        max_tp = tp

    #############################################################################################
    #                                       Media import reactions                              #
    #############################################################################################

    log_and_message("Gap filling from MEDIA", stderr=verbose)
    media_reactions = set()
    for media in growth_media:
        media_reactions.update(PyFBA.gapfill.suggest_from_media(modeldata, reactions, media, verbose))
    added_reactions.append(("media", media_reactions))
    reactions.update(media_reactions)

    for r in media_reactions:
        if r not in reaction_source:
            reaction_source[r] = 'media_reactions'

    tp = measure_accuracy('Media reactions',  growth_media, no_growth_media, reactions, added_reactions, biomass_eqtn,
                          min_growth_conditions, reaction_source, output, verbose)
    if tp > max_tp:
        best_reactions = copy.deepcopy(reactions)
        max_tp = tp

    #############################################################################################
    #                                        Other genomes and organisms                        #
    #############################################################################################

    log_and_message("Gap filling from CLOSE GENOMES", stderr=verbose)
    if close:
        for close_genome in close:
            # add reactions from roles in close genomes
            close_reactions = PyFBA.gapfill.suggest_from_roles(close_genome, modeldata.reactions, threshold=0,
                                                               verbose=verbose)
            # find the new reactions
            close_reactions.difference_update(reactions)
            added_reactions.append((f"close genome: {close_genome}", close_reactions))
            reactions.update(close_reactions)
            for r in close_reactions:
                if r not in reaction_source:
                    reaction_source[r] = f"close genome: {close_genome}"

            tp = measure_accuracy(f"Close genome: {close_genome}",  growth_media, no_growth_media, reactions,
                                  added_reactions, biomass_eqtn, min_growth_conditions, reaction_source, output, verbose)
            if tp > max_tp:
                best_reactions = copy.deepcopy(reactions)
                max_tp = tp

    #############################################################################################
    #                                        Subsystems                                         #
    #############################################################################################

    log_and_message("Gap filling from SUBSYSTEMS", stderr=verbose)
    subsystem_reactions = PyFBA.gapfill.suggest_reactions_from_subsystems(modeldata.reactions, reactions,
                                                                          organism_type=genome_type,
                                                                          threshold=0.5, verbose=verbose)
    added_reactions.append(("subsystems", subsystem_reactions))
    reactions.update(subsystem_reactions)
    for r in subsystem_reactions:
        if r not in reaction_source:
            reaction_source[r] = 'subsystem_reactions'

    tp = measure_accuracy('Subsystems',  growth_media, no_growth_media, reactions, added_reactions, biomass_eqtn,
                          min_growth_conditions, reaction_source, output, verbose)
    if tp > max_tp:
        best_reactions = copy.deepcopy(reactions)
        max_tp = tp

    #############################################################################################
    #                                        Orphan compounds                                   #
    #############################################################################################

    log_and_message("Gap filling from ORPHAN COMPOUNDS", stderr=verbose)
    orphan_compounds = PyFBA.gapfill.suggest_by_compound(modeldata, reactions, 1)
    added_reactions.append(("orphans", orphan_compounds))
    reactions.update(orphan_compounds)
    for r in orphan_compounds:
        if r not in reaction_source:
            reaction_source[r] = 'orphan_compounds'

    tp = measure_accuracy('Orphan compounds',  growth_media, no_growth_media, reactions, added_reactions, biomass_eqtn,
                          min_growth_conditions, reaction_source, output, verbose)
    if tp > max_tp:
        best_reactions = copy.deepcopy(reactions)
        max_tp = tp

    return max_tp, best_reactions


def gapfill_multiple_media():
    orgtypes = ['gramnegative', 'grampositive', 'microbial', 'mycobacteria', 'plant']
    parser = argparse.ArgumentParser(description='Import a list of reactions and then iterate through our gapfilling'
                                                 ' steps to see when we get growth. You can specify multiple --positive'
                                                 ' & --negative media conditions')
    parser.add_argument('-r', '--reactions', help='reactions file', required=True)
    parser.add_argument('-o', '--output', help='file to save new reaction list to', required=True)
    parser.add_argument('-p', '--positive', help='media file(s) on which the organism can grow',
                        required=True, action='append')
    parser.add_argument('-n', '--negative', help='media file(s) on which the organism can NOT grow',
                        default=[], action='append')
    parser.add_argument('-f', '--fraction', help='fraction of growth conditions on which we want growth for success',
                        default=0.8, type=float)
    parser.add_argument('-c', '--close_genomes', help='close genomes reactions file. Multiple files are allowed',
                        action='append')
    parser.add_argument('-t', '--type', default='gramnegative',
                        help=f'organism type for the model (currently allowed are {orgtypes}). Default=gramnegative')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args(sys.argv[2:])

    if not os.path.exists(args.reactions):
        sys.stderr.write(f"FATAL: {args.reactions} does not exist. Please check your files\n")
        sys.exit(1)

    if args.type not in orgtypes:
        sys.exit("Sorry, {} is not a valid organism type".format(args.type))

    log_and_message(f"Running PyFBA with the parameters: {sys.argv}\n", quiet=True)

    # read the enzyme data
    # compounds, reactions, enzymes = PyFBA.parse.model_seed.compounds_reactions_enzymes(args.type)
    global modeldata
    modeldata = PyFBA.parse.model_seed.parse_model_seed_data(args.type, verbose=args.verbose)
    reactions = read_reactions(args.reactions, args.verbose)
    log_and_message(f"Found {len(reactions)} reactions", stderr=args.verbose)

    max_tp, best_reactions = multiple_gapfill(reactions, args.positive, args.negative, args.fraction,
                                              args.close_genomes, args.type, args.output, args.verbose)

    msg = f"Sorry, we added {len(best_reactions)} reactions, but we can never get more than {max_tp} true positives\n"
    msg += f"We have written those reactions to {args.output}"
    with open(args.output, 'w') as out:
        for r in best_reactions:
            out.write(f"{r}\n")




def gapfill_two_media():
    orgtypes = ['gramnegative', 'grampositive', 'microbial', 'mycobacteria', 'plant']
    parser = argparse.ArgumentParser(description='Import a list of reactions and then iterate through our gapfilling'
                                                 ' steps to see when we get growth. You can specify a single --growth'
                                                 ' & --nogrowth media conditions')
    parser.add_argument('-r', '--reactions', help='reactions file', required=True)
    parser.add_argument('-o', '--output', help='file to save new reaction list to', required=True)
    parser.add_argument('-g', '--growth', help='media file on which the organism can grow',
                        required=True)
    parser.add_argument('-n', '--nogrowth', help='media file on which the organism can NOT grow',
                        default=[])
    parser.add_argument('-c', '--close_genomes', help='close genomes reactions file. Multiple files are allowed',
                        action='append')
    parser.add_argument('-t', '--type', default='gramnegative',
                        help=f'organism type for the model (currently allowed are {orgtypes}). Default=gramnegative')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args(sys.argv[2:])

    if not os.path.exists(args.reactions):
        sys.stderr.write(f"FATAL: {args.reactions} does not exist. Please check your files\n")
        sys.exit(1)

    if args.type not in orgtypes:
        sys.exit("Sorry, {} is not a valid organism type".format(args.type))

    log_and_message(f"Running PyFBA with the parameters: {sys.argv}\n", quiet=True)

    # read the enzyme data
    # compounds, reactions, enzymes = PyFBA.parse.model_seed.compounds_reactions_enzymes(args.type)
    global modeldata
    modeldata = PyFBA.parse.model_seed.parse_model_seed_data(args.type, verbose=args.verbose)
    reactions = read_reactions(args.reactions, args.verbose)
    log_and_message(f"Found {len(reactions)} reactions", stderr=args.verbose)

    growth_media = PyFBA.parse.read_media.find_media_file(args.growth, modeldata=modeldata, verbose=args.verbose)
    no_growth_media = PyFBA.parse.read_media.find_media_file(args.nogrowth, modeldata=modeldata, verbose=args.verbose)


    PyFBA.gapfill.gapfill_two_media(modeldata, reactions, growth_media, no_growth_media,
                                                  args.close_genomes, args.type, args.output, verbose=args.verbose)
