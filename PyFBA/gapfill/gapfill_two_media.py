"""
Gapfill given two media - something where we grow and something where we don't grow.

We gapfill like normal until we grow, and then bisect the growth, but we only
allow it to proceed if we don't grow on the negative media.
"""
import copy
import os
import sys
import argparse

import PyFBA
from PyFBA import log_and_message


def gapfill_two_media(model_data, reactions, growth_media, no_growth_media, close, genome_type, output, verbose=False):
    """
    Run multiple gap filling operations and try to resove positive/negative growth
    :param model_data: The SEED Model Data
    :param output: Output file name
    :type output: str
    :param reactions: A list of reactions we think are in the genome
    :type reactions: set[str]
    :param growth_media: The media object where we grow
    :type growth_media: set
    :param no_growth_media: The media object where we DO NOT grow
    :type no_growth_media: set
    :param close: A list of close genomes used for gapfilling
    :type close: list[str]
    :param genome_type: The genome type (gram +ve, -ve, etc)
    :type genome_type: str
    :param output: the name of the file to write the reactions to
    :type output: str
    :param verbose: more output
    :type verbose: bool
    :return: The set of reactions that meets the conditions
    :rtype: set[str]
    """

    biomass_eqtn = PyFBA.metabolism.biomass.biomass_equation(genome_type)

    # initial check. Do we NOT grow in no growth media?
    status, value, growth = PyFBA.fba.run_fba(model_data, reactions, no_growth_media, biomass_eqtn, verbose=verbose)
    if growth:
        log_and_message(f"There was growth on the NO GROWTH media.", stderr=True, loglevel="WARN")
        sys.exit(-1)

    # initial check. Do we NOT grow in growth media?
    status, value, growth = PyFBA.fba.run_fba(model_data, reactions, growth_media, biomass_eqtn, verbose=verbose)
    if growth:
        log_and_message(f'We already grow on the GROWTH media. There is nothing to do!', stderr=True, loglevel='WARN')
        sys.exit(0)

    # gap fill this until we get growth
    gfreactions = PyFBA.gapfill.gapfill(reactions, model_data, growth_media, biomass_eqtn,
                                                        close, genome_type, verbose=verbose)
    # now we need to check again!
    status, value, nogrowth = PyFBA.fba.run_fba(model_data, reactions, no_growth_media, biomass_eqtn, verbose=verbose)
    status, value, growth = PyFBA.fba.run_fba(model_data, reactions, growth_media, biomass_eqtn, verbose=verbose)

    count = 0
    exclude_reactions = set()
    while growth and nogrowth and count < 5:
        count += 1
        # if we are here we are going to iterate through one at a time and remove any reaction that
        # allows us to grow on the nogrowth media, and then gapfill again
        allreactions = list(gfreactions.keys())
        for idx in range(len(allreactions)):
            testr = allreactions[:idx] + allreactions[idx+1:]
            status, value, growth = PyFBA.fba.run_fba(model_data, set(testr), no_growth_media, biomass_eqtn,
                                                      verbose=verbose)
            if not growth:
                # this reaction allows us to grow in this condition, so exclude it!
                log_and_message(f"Reaction {allreactions[idx]} allows us to grow where we shouldn't. Skipped", stderr=verbose)
                exclude_reactions.add(allreactions[idx])

        gfreactions = PyFBA.gapfill.gapfill(reactions, model_data, growth_media, biomass_eqtn,
                                            close, genome_type, r2exclude=exclude_reactions, verbose=verbose)
        # now we need to check again!
        status, value, nogrowth = PyFBA.fba.run_fba(model_data, reactions, no_growth_media, biomass_eqtn,
                                                    verbose=verbose)
        status, value, growth = PyFBA.fba.run_fba(model_data, reactions, growth_media, biomass_eqtn, verbose=verbose)

    if growth and not nogrowth:
        log_and_message(f'Complete. We are able to grow on the growth media and not on the no growth media', stderr=verbose)
        log_and_message(f'Reaction list has been written to {output}', stderr=verbose)
        with open(output, 'w') as out:
            for r in gfreactions:
                out.write(f"{r}\t{gfreactions[r]}\n")
        exit(0)
    else:
        log_and_message(f"ERROR: We were never able to get growth and no growth correct. Gave up trying after {count} iterations!")
        log_and_message(f"Nonetheless we wrote the reactions to {output}")
        with open(output, 'w') as out:
            out.write("# PYFBA Warning: We could not resolve whether growth/nogrowth is correct for these reactions")
            for r in gfreactions:
                out.write(f"{r}\t{gfreactions[r]}\n")
        exit(-1)