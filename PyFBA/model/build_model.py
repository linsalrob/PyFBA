from __future__ import print_function
import sys
import PyFBA


def roles_to_model(rolesFile, id, name, orgtype="gramnegative", verbose=False):
    """
    Read in the 'assigned_functions' file from RAST and create a model.

    :param rolesFile: File path to assigned functions RAST file
    :type rolesFile: str
    :param id: Model ID
    :type id: str
    :param name: Model name
    :type name: str
    :param orgtype: Organism type
    :type orgtype: str
    :param verbose: Verbose output
    :type verbose: bool
    :return: The generated model object
    :rtype: Model
    """

    # Load ModelSEED database
    compounds, reactions, enzymes = \
            PyFBA.parse.model_seed.compounds_reactions_enzymes(orgtype)

    # Read in assigned functions file to build set of roles
    assigned_functions = PyFBA.parse.read_assigned_functions(rolesFile)
    roles = set()
    for rs in assigned_functions.values():
        roles.update(rs)

    # Obtain reactions for each role
    # Key is role, value is set of reaction ids
    model_reactions = PyFBA.filters.roles_to_reactions(roles)

    # Create model object
    model = PyFBA.model.Model(id, name)
    for role, rxnIDs in model_reactions.items():
        for rxnID in rxnIDs:
            if rxnID in reactions:
                model.add_reactions({reactions[rxnID]})
                model.add_roles({role: {rxnID}})
            elif verbose:
                print("Reaction ID '{}' for role '{}'".format(rxnID, role),
                      "is not in our reactions list. Skipped.",
                      file=sys.stderr)

    # Set biomass equation based on organism type
    biomass_eqn = PyFBA.metabolism.biomass_equation(orgtype)
    model.set_biomass_reaction(biomass_eqn)
    return model
