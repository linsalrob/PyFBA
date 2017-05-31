from __future__ import print_function
import sys
import os.path
import PyFBA


def common_compounds():
    """
    Return set of highly common compounds
    """
    cc = ["H+", "H2O", "ATP", "ADP","Phosphate",
          "NAD", "NADH", "NADP", "NADPH", "PPi",
          "CoA", "CO2"]
    cpds = set()
    for c in cc:
        cpds.add(PyFBA.metabolism.Compound(c, "c"))
        cpds.add(PyFBA.metabolism.Compound(c, "e"))
    return cpds


def rxn_to_rxn_network(model, filter_threshold=100):
    """
    Build and return a simple reaction-to-reaction interaction mapping.
    Mapping will be in a Dictionary of:
        rxnID -> set(list of rxnIDs)
        key: String
        value: set(String)

    :param model: Model object to obtain interaction mapping from
    :type model: Model
    :param filter_threshold: compounds with this amount will be filtered out
                             when finding reaction connections. Set to -1 if
                             no threshold is to be used
    :type filter_threshold: int
    :rtype: dict
    """
    network = {}  # Will hold full network interactions

    # Get list of most compounds to ignore
    skip = common_compounds()

    rxnIDset = set(model.reactions.keys())
    # Iterate through reactions in model
    for rID, rxn in model.reactions.items():
        network[rID] = set()
        # Get compounds associated with this reaction
        cpds = rxn.left_compounds.union(rxn.right_compounds)
        # Remove compounds from skip set
        cpds -= skip

        # Only keep compounds that are below filter threshold
        filtered_cpds = set()
        for c in cpds:
            if filter_threshold != -1 and len(c.reactions) < filter_threshold:
                filtered_cpds.add(c)
        cpds = filtered_cpds

        # For each compound, get list of reactions
        for c in cpds:
            cpd_rxns = c.reactions
            # Grab only reactions that are in the current model
            cpd_rxns &= rxnIDset
            # Make sure its own reaction ID is not included here
            cpd_rxns -= set([rID])
            # Add reactions to network
            network[rID].update(cpd_rxns)

    return network


def number_of_nodes(network):
    """
    Return number of nodes in a network

    param network: Network generated from rxn_to_rxn_network()
    type network: dict
    rtype: int
    """
    return len(network)


def clustering_coeff(network):
    pass
