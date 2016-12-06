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


def stdout_fba(model, media_file, biomass_reaction=None):
    """
    Run FBA on model and output results in tab-delimited format.

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
    print("reaction\tflux\tfunction")
    for r, flux in fluxes.items():
        if r not in mReactions or mReactions[r] == "":
            rolecolumn = "None"
        else:
            rolecolumn = ";".join(mReactions[r])
        print(r, flux, rolecolumn)


def stdout_model(model):
    """
    Output model reaction, function, and gap-fill information.

    :param model: Model object
    :type model: Model
    """
    # Get mapping from reaction IDs to roles
    mReactions = {r: [] for r in model.reactions.keys()}
    for role, rxns in model.roles.items():
        for r in rxns:
            mReactions[r].append(role)
    # Print header
    print("reaction\tfunction\tequation\tgapfilled")

    # Print reactions from model
    for r, roles in mReactions.items():
        eqn = model.reactions[r].equation
        rolecolumn = ";".join(roles)
        print(r, rolecolumn, eqn, "no", sep="\t")

    # Print reactions from gap-filling
    for r, roles in model.gf_reactions.items():
        eqn = model.reactions[r].equation
        rolecolumn = ";".join(roles)
        print(r, rolecolumn, eqn, "yes", sep="\t")
