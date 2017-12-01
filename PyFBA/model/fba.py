from __future__ import print_function
import sys
import os.path
import PyFBA


def model_reaction_fluxes(model, media_file, biomass_reaction=None):
    """
    Run FBA on model and return dictionary of reaction ID and flux.

    :param model: Model object to obtain fluxes from
    :type model: PyFBA.model.Model
    :param media_file: Media filepath
    :type media_file: str
    :param biomass_reaction: Given biomass Reaction object
    :type biomass_reaction: PyFBA.metabolism.Reaction
    :rtype: dict
    """
    status, value, growth = model.run_fba(media_file, biomass_reaction)
    if not growth:
        print("Warning: model did not grow on given media", file=sys.stderr)
    return PyFBA.fba.reaction_fluxes()


def output_fba(f, model, media_file, biomass_reaction=None):
    """
    Run FBA on model and output results in tab-delimited format.

    :param f: File object to print to
    :type f: file
    :param model: Model object to obtain fluxes from
    :type model: PyFBA.model.Model
    :param media_file: Media filepath
    :type media_file: str
    :param biomass_reaction: Given biomass Reaction object
    :type biomass_reaction: PyFBA.metabolism.Reaction
    """
    # Get mapping from reaction IDs to roles
    mReactions = {r: [] for r in model.reactions.keys()}
    for role, rxns in model.roles.items():
        for r in rxns:
            mReactions[r].append(role)

    # Run FBA and get fluxes
    fluxes = model_reaction_fluxes(model, media_file, biomass_reaction)

    # Print header
    f.write("reaction\tflux\tfunction\n")
    for r, flux in fluxes.items():
        if r not in mReactions or mReactions[r] == "":
            rolecolumn = "None"
        else:
            rolecolumn = ";".join(mReactions[r])
        f.write("\t".join([r, str(flux), rolecolumn]))
        f.write("\n")


def output_fba_with_subsystem(f, model, media_file, biomass_reaction=None):
    """
    Run FBA on model and output results and subsystem info in tab-delimited format.

    :param f: File object to print to
    :type f: file
    :param model: Model object to obtain fluxes from
    :type model: PyFBA.model.Model
    :param media_file: Media filepath
    :type media_file: str
    :param biomass_reaction: Given biomass Reaction object
    :type biomass_reaction: PyFBA.metabolism.Reaction
    """
    # Get mapping from reaction IDs to roles
    mReactions = {r: [] for r in model.reactions.keys()}
    for role, rxns in model.roles.items():
        for r in rxns:
            mReactions[r].append(role)

    # Load subsystem info
    ss_data = {}
    with open(os.path.join(os.path.dirname(__file__), "..", "util", "full_roles_ss.tsv")) as fin:
        for l in fin:
            func, cat, subcat, ss = l.rstrip("\n").split("\t")
            if func not in ss_data:
                ss_data[func] = set()
            cat = cat if cat != "" else "Unknown"
            subcat = subcat if subcat != "" else "Unknown"
            ss = ss if ss != "" else "Unknown"
            ss_data[func].add((cat, subcat, ss))

    # Run FBA and get fluxes
    fluxes = model_reaction_fluxes(model, media_file, biomass_reaction)

    # Print header
    f.write("reaction\tflux\tfunction\tsubsystem\tsubcategory\tcategory\n")
    # Iterate through reactions and their fluxes
    for r, flux in fluxes.items():
        try:
            # Iterate through the roles for each reaction
            for role in mReactions[r]:
                try:
                    for info in ss_data[role]:
                        cat, subcat, ss = info
                        f.write("\t".join([r, str(flux), role, ss, subcat, cat]))
                        f.write("\n")
                except KeyError:
                    # This means we have no subsystem information for the role
                    cat = subcat = ss = "Unknown"
                    f.write("\t".join([r, str(flux), role, ss, subcat, cat]))
                    f.write("\n")
        except KeyError:
            # This means we have no role information for the reaction
            role = ss = subcat = cat = "Unknown"
            f.write("\t".join([r, str(flux), role, ss, subcat, cat]))
            f.write("\n")


def compute_compound_counts(model, media_file, biomass_reaction=None,
                            verbose=False):
    """
    Using FBA, calculate overall compound counts after summing reaction fluxes

    :param model: Model object
    :type model: PyFBA.model.Model
    :param media_file: Media filepath
    :type media_file: str
    :param biomass_reaction: Given biomass Reaction object
    :type biomass_reaction: PyFBA.metabolism.Reaction
    :return: Compound counts
    :rtype: dict
    """
    # Run FBA and obtain reaction fluxes
    fluxes = model_reaction_fluxes(model, media_file, biomass_reaction)
    if verbose:
        print("FBA completed", file=sys.stderr)

    # Load ModelSEED database
    compounds, reactions, enzymes = \
        PyFBA.parse.compounds_reactions_enzymes(model.organism_type)

    cpd_counts = {"e": {}, "c": {}}
    # Iterate through fluxes to calculate compound counts
    for rxnID, rflux in fluxes.items():
        # First check for biomass reactions
        if rxnID == "BIOMASS_EQN":
            cpd_counts["c"]["Biomass"] = rflux
            continue

        # Biomass has its own secretion reaction, we will skip it
        elif rxnID == "UPTAKE_SECRETION_REACTION cpd11416":
            continue

        # Uptake and secretion reactions will be included here but not
        # present in the ModelSEED database
        # The compound ModelSEED ID is located at the end of the name
        elif rxnID.startswith("UPTAKE_SECRETION_REACTION"):
            cpdID = rxnID.replace("UPTAKE_SECRETION_REACTION ", "")

            # Find the compound in the model
            cpd = None
            for c in model.compounds:
                if c.model_seed_id == cpdID and c.location == "e":
                    cpd = c
                    break

            # Check if compound was found
            if cpd is None and rflux != 0 and verbose:
                print("Compound {} was not found in model but was found in an "
                      "uptake and secretion reaction with flux".format(cpdID),
                      file=sys.stderr)
                print("Reaction: {}\tFlux:{}".format(rxnID, rflux),
                      file=sys.stderr)
                continue

            # Add flux to the cpd_count
            cpd_name = str(cpd)
            if cpd_name not in cpd_counts["e"]:
                cpd_counts["e"][cpd_name] = 0

            # The right side of these reactions are always the boundary
            # sink and the left side are always the external compound
            # Negate flux since a negative flux of this reaction means
            # the compound is being produced and not consumed
            cpd_counts["e"][cpd_name] += -rflux
            continue

        # Get left compounds
        for cpd in reactions[rxnID].left_compounds:
            loc = cpd.location
            cpd_name = str(cpd)
            if cpd_name not in cpd_counts[loc]:
                cpd_counts[loc][cpd_name] = 0

            # Add flux to cpd_count
            # Check direction because a leftward reaction means left
            # compounds are products and not reactants
            if reactions[rxnID].direction == "<":
                cpd_counts[loc][cpd_name] += rflux
            else:
                cpd_counts[loc][cpd_name] += -rflux

        # Get right compounds
        for cpd in reactions[rxnID].right_compounds:
            cpd_name = str(cpd)
            loc = cpd.location
            if cpd_name not in cpd_counts[loc]:
                cpd_counts[loc][cpd_name] = 0

            # Add flux to cpd_count
            # Check direction because a leftward reaction means right
            # compounds are reactants and not products
            if reactions[rxnID].direction == "<":
                cpd_counts[loc][cpd_name] += -rflux
            else:
                cpd_counts[loc][cpd_name] += rflux

    return cpd_counts
