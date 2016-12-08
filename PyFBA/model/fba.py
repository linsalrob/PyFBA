from __future__ import print_function
import sys
import PyFBA


def model_reaction_fluxes(model, media_file, biomass_reaction=None):
    """
    Run FBA on model and return dictionary of reaction ID and flux.

    :param model: Model object to obtain fluxes from
    :type model: Model
    :param media_file: Media filepath
    :type media_file: str
    :param biomass_reaction: Given biomass Reaction object
    :type biomass_reaction: Reaction
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
    :type model: Model
    :param media_file: Media filepath
    :type media_file: str
    :param biomass_reaction: Given biomass Reaction object
    :type biomass_reaction: Reaction
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
        f.write("\t".join([r, flux, rolecolumn]))
        f.write("\n")
