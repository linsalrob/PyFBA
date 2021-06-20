"""
Given a set of roles (e.g. from a genome annotation) can we gap fill those?
Largely based on From_functional_roles_to_gap-filling
"""

import os
import sys
import PyFBA
import argparse
import copy
from PyFBA import log_and_message


def run_eqn(why, md, r2r, med, bme, verbose=False):
    """
    Run the fba
    :param why: why are we doing this
    :param md: modeldata
    :param r2r: reactions to run
    :param med: media object
    :param bme: biomass equation
    :param verbose: more output
    :type verbose: bool
    :return: (value, growth)
    """

    status, value, growth = PyFBA.fba.run_fba(md, r2r, med, bme)
    log_and_message(f"FBA run {why} has a biomass flux value of {value} --> Growth: {growth}", stderr=verbose)
    return value, growth


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
        new_essential = PyFBA.gapfill.minimize_additional_reactions(ori, new, modeldata, media, biomass_equation,
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


def roles_to_reactions_to_run(roles, orgtype='gramnegative', verbose=False):
    roles_to_reactions = PyFBA.filters.roles_to_reactions(roles, organism_type=orgtype, verbose=verbose)
    reactions_to_run = set()
    for role in roles_to_reactions:
        reactions_to_run.update(roles_to_reactions[role])
    log_and_message(f"There are {len(reactions_to_run)} unique reactions associated with this genome", stderr=verbose)
    return reactions_to_run


def read_media(mediafile, modeldata, verbose=False):
    """
    Read the media file and return a set of compounds
    :param modeldata: the modeldata object
    :type modeldata: PyFBA.model_seed.ModelData
    :param mediafile: the media file to read
    :param verbose: more output
    :type verbose: bool
    :return: a set of media compounds
    :rtype: Set[PyFBA.metabolism.Compound]
    """

    if mediafile in PyFBA.parse.media_files():
        log_and_message(f"parsing media directly from {mediafile}", stderr=verbose)
        # pyfba media already corrects the names, so we can  just return it.
        return PyFBA.parse.pyfba_media(mediafile, modeldata)
    elif os.path.exists(mediafile):
        log_and_message(f"parsing media file {mediafile}", stderr=verbose)
        media = PyFBA.parse.read_media_file(mediafile)
    elif 'PYFBA_MEDIA_DIR' in os.environ and os.path.exists(os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediafile)):
        log_and_message(f"parsing media file {os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediafile)}", stderr=verbose)
        media = PyFBA.parse.read_media_file(os.path.join(os.environ['PYFBA_MEDIA_DIR'], mediafile))
    else:
        log_and_message(f"Can't figure out how to parse media from {mediafile}", stderr=True, loglevel="CRITICAL")
        sys.exit(-1)
    return PyFBA.parse.correct_media_names(media, modeldata.compounds)


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


def run_gapfill_from_roles(roles, reactions_to_run, modeldata, media, orgtype='gramnegative', close_orgs=None,
                           close_genera=None, verbose=False):
    """
    gapfill growth from a set of roles in the genome
    :param close_genera: the list of roles in close genera
    :param close_orgs: the list of roles in close organisms
    :param roles: The set of roles in this genome
    :type roles: set[str[
    :param reactions_to_run: The reactions to run
    :type reactions_to_run: set[str]
    :param modeldata: the modeldata object
    :type modeldata: PyFBA.model_seed.ModelData
    :param media:  a set of media compounds
    :type media: Set[PyFBA.metabolism.Compound]
    :param orgtype: the organism type for the model
    :type orgtype: str
    :param verbose: more output
    :type verbose: bool
    :return: a dict of the reactions and what step they were added at
    :rtype: dict[str, str]
    """

    tempset = set()
    for r in reactions_to_run:
        if r in modeldata.reactions:
            tempset.add(r)
        else:
            log_and_message(f"Reaction ID {r} is not in our reactions list. Skipped", stderr=verbose)
    reactions_to_run = tempset

    biomass_equation = PyFBA.metabolism.biomass_equation(orgtype)

    run_eqn("Initial", modeldata, reactions_to_run, media, biomass_equation, verbose=verbose)

    added_reactions = []
    original_reactions_to_run = copy.deepcopy(reactions_to_run)

    #############################################################################################
    #                                       Gapfilling                                          #
    #                                                                                           #
    #  We do this in the order:                                                                 #
    #     essential reactions: because you need to have these, but it is stronger evidence if   #
    #              your friends have it too!                                                    #
    #     media: because you should be importing everything in the media                        #
    #     linked_reactions: because they make sense!                                            #
    #     closely related organisms: because you should have roles your friends have            #
    #     subsystems: to complete things you already have                                       #
    #     orphans: to make sure everything is produced/consumed                                 #
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
    reactions_to_run = update_r2r(reactions_to_run, essential_reactions, "ESSENTIAL REACTIONS")
    value, growth = run_eqn("Initial", modeldata, reactions_to_run, media, biomass_equation, verbose=verbose)

    if growth:
        return minimize_reactions(original_reactions_to_run, added_reactions, modeldata, media, biomass_equation,
                                  verbose=verbose)

    #############################################################################################
    #                                       LINKED REACTIONS                                    #
    #############################################################################################

    log_and_message("Gap filling from Linked Reactions", stderr=verbose)
    linked_reactions = PyFBA.gapfill.suggest_linked_reactions(modeldata, reactions_to_run)
    for r in linked_reactions:
        modeldata.reactions[r].reset_bounds()
    added_reactions.append(("linked_reactions", linked_reactions))
    reactions_to_run = update_r2r(reactions_to_run, linked_reactions, "LINKED REACTIONS")
    value, growth = run_eqn("Initial", modeldata, reactions_to_run, media, biomass_equation, verbose=verbose)

    if growth:
        return minimize_reactions(original_reactions_to_run, added_reactions, modeldata, media, biomass_equation,
                                  verbose=verbose)

    #############################################################################################
    #                                       EC NUMBERS                                          #
    #############################################################################################

    log_and_message("Gap filling from limited EC numbers", stderr=verbose)
    ecnos = PyFBA.gapfill.suggest_reactions_using_ec(roles, modeldata, reactions_to_run, verbose=verbose)
    for r in ecnos:
        modeldata.reactions[r].reset_bounds()
    added_reactions.append(("ec_numbers_brief", ecnos))
    reactions_to_run = update_r2r(reactions_to_run, ecnos, "EC Numbers")
    value, growth = run_eqn("Initial", modeldata, reactions_to_run, media, biomass_equation, verbose=verbose)

    if growth:
        return minimize_reactions(original_reactions_to_run, added_reactions, modeldata, media, biomass_equation,
                                  verbose=verbose)

    #############################################################################################
    #                                       Media import reactions                              #
    #############################################################################################

    log_and_message("Gap filling from MEDIA", stderr=verbose)
    media_reactions = PyFBA.gapfill.suggest_from_media(modeldata, reactions_to_run, media, verbose=verbose)
    added_reactions.append(("media", media_reactions))
    reactions_to_run = update_r2r(reactions_to_run, media_reactions, "MEDIA REACTIONS")
    value, growth = run_eqn("Initial", modeldata, reactions_to_run, media, biomass_equation, verbose=verbose)

    if growth:
        return minimize_reactions(original_reactions_to_run, added_reactions, modeldata, media, biomass_equation,
                                  verbose=verbose)

    #############################################################################################
    #                                        Other genomes and organisms                        #
    #############################################################################################

    log_and_message("Gap filling from CLOSE GENOMES", stderr=verbose)
    if close_orgs:
        # add reactions from roles in close genomes
        close_reactions = PyFBA.gapfill.suggest_from_roles(close_orgs, modeldata.reactions, threshold=0,
                                                           verbose=verbose)
        close_reactions.difference_update(reactions_to_run)
        added_reactions.append(("close genomes ", close_reactions))
        reactions_to_run = update_r2r(reactions_to_run, close_reactions, "CLOSE ORGANISMS")
        value, growth = run_eqn("Initial", modeldata, reactions_to_run, media, biomass_equation, verbose=verbose)

        if growth:
            return minimize_reactions(original_reactions_to_run, added_reactions, modeldata, media, biomass_equation,
                                      verbose=verbose)

    if close_genera:
        # add reactions from roles in similar genera
        genus_reactions = PyFBA.gapfill.suggest_from_roles(close_genera, modeldata.reactions, threshold=0,
                                                           verbose=verbose)
        genus_reactions.difference_update(reactions_to_run)
        added_reactions.append(("other genera", genus_reactions))
        reactions_to_run = update_r2r(reactions_to_run, genus_reactions, "CLOSE GENERA")
        value, growth = run_eqn("Initial", modeldata, reactions_to_run, media, biomass_equation, verbose=verbose)

        if growth:
            return minimize_reactions(original_reactions_to_run, added_reactions, modeldata, media, biomass_equation,
                                      verbose=verbose)

    #############################################################################################
    #                                        Subsystems                                         #
    #############################################################################################

    log_and_message("Gap filling from SUBSYSTEMS", stderr=verbose)
    subsystem_reactions = PyFBA.gapfill.suggest_reactions_from_subsystems(modeldata.reactions, reactions_to_run,
                                                                          organism_type=orgtype, threshold=0.5,
                                                                          verbose=verbose)
    added_reactions.append(("subsystems", subsystem_reactions))
    reactions_to_run = update_r2r(reactions_to_run, subsystem_reactions, "SUBSYSTEMS")
    value, growth = run_eqn("Initial", modeldata, reactions_to_run, media, biomass_equation, verbose=verbose)

    if growth:
        return minimize_reactions(original_reactions_to_run, added_reactions, modeldata, media, biomass_equation,
                                  verbose=verbose)
    #############################################################################################
    #                                        Orphan compounds                                   #
    #############################################################################################

    log_and_message("Gap filling from ORPHANS", stderr=verbose)
    orphan_compounds = PyFBA.gapfill.suggest_by_compound(modeldata, reactions_to_run, 1)
    added_reactions.append(("orphans", orphan_compounds))
    reactions_to_run = update_r2r(reactions_to_run, orphan_compounds, "ORPHANS")
    value, growth = run_eqn("Initial", modeldata, reactions_to_run, media, biomass_equation, verbose=verbose)

    if growth:
        return minimize_reactions(original_reactions_to_run, added_reactions, modeldata, media, biomass_equation,
                                  verbose=verbose)

    # ## Revisit EC Numbers
    #
    # When we added the EC numbers before, we were a little conservative, only adding those EC numbers that appeared in
    # two or less (by default) reactions. If we get here, lets be aggressive and add any EC number regardless of how
    # many reactions we add. We set the `maxnumrx` variable to 0
    #############################################################################################
    #                                       EC NUMBERS                                          #
    #############################################################################################

    log_and_message("Gap filling from limited EC numbers", stderr=verbose)
    ecnos = PyFBA.gapfill.suggest_reactions_using_ec(roles, modeldata, reactions_to_run, maxnumrx=0, verbose=verbose)
    for r in ecnos:
        modeldata.reactions[r].reset_bounds()
    added_reactions.append(("ec_numbers_full", ecnos))
    reactions_to_run = update_r2r(reactions_to_run, ecnos, "EC Numbers")
    value, growth = run_eqn("Initial", modeldata, reactions_to_run, media, biomass_equation, verbose=verbose)

    if growth:
        return minimize_reactions(original_reactions_to_run, added_reactions, modeldata, media, biomass_equation,
                                  verbose=verbose)

    # We revist linked reactions once more, because now we have many more reactions in our set to run!

    #############################################################################################
    #                                       LINKED REACTIONS                                    #
    #############################################################################################

    log_and_message("Gap filling from Linked Reactions", stderr=verbose)
    linked_reactions = PyFBA.gapfill.suggest_linked_reactions(modeldata, reactions_to_run)
    for r in linked_reactions:
        modeldata.reactions[r].reset_bounds()
    added_reactions.append(("linked_reactions_full", linked_reactions))
    reactions_to_run = update_r2r(reactions_to_run, linked_reactions, "LINKED REACTIONS")
    value, growth = run_eqn("Initial", modeldata, reactions_to_run, media, biomass_equation, verbose=verbose)

    if growth:
        return minimize_reactions(original_reactions_to_run, added_reactions, modeldata, media, biomass_equation,
                                  verbose=verbose)

    log_and_message(f"FATAL: After compiling {len(reactions_to_run)} reactions, we still could not get growth",
                    stderr=True, loglevel='CRITICAL')
    return set()


def gapfill_from_roles():
    """
    Parse the arguments and start the gapfilling.
    """

    orgtypes = ['gramnegative', 'grampositive', 'microbial', 'mycobacteria', 'plant']
    parser = argparse.ArgumentParser(description='Run Flux Balance Analysis on a set of gapfilled functional roles')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-r', '--roles', help='A list of functional roles in this genome, one per line')
    group.add_argument('-a', '--assigned_functions', help='RAST assigned functions (tab separated PEG/Functional Role)')
    group.add_argument('-f', '--features', help='PATRIC features.txt file (with 5 columns)')
    parser.add_argument('-o', '--output', help='file to save new reaction list to', required=True)
    parser.add_argument('-m', '--media', help='media name', required=True)
    parser.add_argument('-t', '--type', default='gramnegative',
                        help=f'organism type for the model (currently allowed are {orgtypes}). Default=gramnegative')
    parser.add_argument('-c', '--close', help='a file with roles from close organisms')
    parser.add_argument('-g', '--genera', help='a file with roles from similar genera')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args(sys.argv[2:])

    log_and_message(f"Running PyFBA with the parameters: {sys.argv}\n", quiet=True)

    model_data = PyFBA.parse.model_seed.parse_model_seed_data(args.type)

    if args.roles:
        if not os.path.exists(args.roles):
            sys.stderr.write(f"FATAL: {args.roles} does not exist. Please check your files\n")
            sys.exit(1)
        log_and_message(f"Getting the roles from {args.roles}", stderr=args.verbose)
        roles = PyFBA.parse.read_functional_roles(args.roles, args.verbose)
    elif args.assigned_functions:
        if not os.path.exists(args.assigned_functions):
            sys.stderr.write(f"FATAL: {args.assigned_functions} does not exist. Please check your files\n")
            sys.exit(1)
        log_and_message(f"Getting the roles from {args.assigned_functions}", stderr=args.verbose)
        roles = PyFBA.parse.assigned_functions_set(args.assigned_functions)
    elif args.features:
        if not os.path.exists(args.features):
            sys.stderr.write(f"FATAL: {args.features} does not exist. Please check your files\n")
            sys.exit(1)
        log_and_message(f"Getting the roles from {args.features}", stderr=args.verbose)
        roles = PyFBA.parse.read_features_file(args.features, args.verbose)
    else:
        sys.stderr.write("FATAL. Either a roles or functions file must be provided")
        sys.exit(1)

    reactions_to_run = roles_to_reactions_to_run(roles, args.type, args.verbose)
    media = read_media(args.media, model_data, args.verbose)

    new_reactions = run_gapfill_from_roles(roles=roles, reactions_to_run=reactions_to_run, modeldata=model_data,
                                           media=media, orgtype=args.type, close_orgs=args.close,
                                           close_genera=args.genera, verbose=args.verbose)
    if new_reactions:
        with open(args.output, 'w') as out:
            for r in new_reactions:
                out.write(f"{r}\t{new_reactions[r]}\n")


if __name__ == "__main__":
    gapfill_from_roles()
