from __future__ import print_function, absolute_import, division
import sys
from itertools import combinations


def network_compounds_sif(network,  filepath, ulimit=None, verbose=False):
    """
    Create a SIF (simple interaction file) format of a network. Network nodes
    are compounds. SIF is compatible with Cytoscape. The filepath provided
    will be overwritten. An upper limit of number of neighbors a compound
    can have is provided.
    Note, the following compounds typically have a high number of connections:
    CO2 (location: c)
    ADP (location: c)
    ACP (location: c)
    NAD (location: c)
    NADH (location: c)
    CoA (location: c)
    NADP (location: c)
    H2O (location: c)
    L-Glutamate (location: c)
    Pyruvate (location: c)
    Phosphate (location: c)
    Glycerol-3-phosphate (location: c)
    O2 (location: c)
    H+ (location: e)
    NH3 (location: c)
    ATP (location: c)
    H+ (location: c)
    PPi (location: c)
    NADPH (location: c)
    AMP (location: c)

    :param network: Network to produce file for
    :type network: PyFBA.network.Network
    :param filepath: Filepath for output
    :type filepath: str
    :param ulimit: For a compound, the highest allowed number of neighbors
                   to be included in the output
    :type ulimit: int
    :param verbose: Verbose output flag
    :type verbose: bool
    :return: None
    """
    fh = open(filepath + ".sif", "w")
    conns = set()  # Record which connections were already passed through
    conns_to_output = set()  # Connections to print out
    counts = dict()

    num_nodes = network.number_of_nodes()
    num_edges = network.number_of_edges()
    if verbose:
        print("Network contains", num_nodes, "nodes and", num_edges, "edges",
              file=sys.stderr)

    # If ulimit is not set, go through edges since it is much faster than
    # iterating through nodes
    if ulimit is None:
        # Iterate through edges
        for i, e in enumerate(network.edges_iter(), start=1):
            if verbose:
                print("Working on edge {} out of {}".format(i, num_edges),
                      end="\r", file=sys.stderr)
            # Skip connections we've already seen
            # Could be due to bidirectional reactions
            if e in conns:
                continue
            cpd1, cpd2 = e
            conns.add(e)
            conns.add((cpd2, cpd1))
            if cpd1 not in counts:
                counts[cpd1] = 0
            if cpd2 not in counts:
                counts[cpd2] = 0
            counts[cpd1] += 1
            counts[cpd2] += 1
            fh.write(cpd1.name + "\tcc\t" + cpd2.name + "\n")
        fh.close()

        with open(filepath + ".sif", "w") as fh:
            for cpd, count in counts.items():
                fh.write(cpd.name + "\t" + count + "\n")

    else:
        # Iterate through edges
        for i, e in enumerate(network.edges_iter(), start=1):
            if verbose:
                print("Working on edge {} out of {}".format(i, num_edges),
                      end="\r", file=sys.stderr)
            # Skip connections we've already seen
            # Could be due to bidirectional reactions
            if e in conns:
                continue
            cpd1, cpd2 = e
            conns.add(e)
            conns.add((cpd2, cpd1))
            conns_to_output.add(e)
            if cpd1 not in counts:
                counts[cpd1] = 0
            if cpd2 not in counts:
                counts[cpd2] = 0
            counts[cpd1] += 1
            counts[cpd2] += 1

        for c in conns_to_output:
            cpd1, cpd2 = c
            if counts[cpd1] > ulimit or counts[cpd2] > ulimit:
                continue
            fh.write(cpd1.name + "\tcc\t" + cpd2.name + "\n")
        fh.close()

        with open(filepath + ".sizes", "w") as fh,\
                open(filepath + ".skip", "w") as fhskip:
            for cpd, count in counts.items():
                if count > ulimit:
                    fhskip.write(cpd.name + "\t" + str(count) + "\n")
                    continue
                fh.write(cpd.name + "\t" + str(count) + "\n")
    if verbose:
        print("", file=sys.stderr)


def network_reactions_sif(network,  filepath):
    """
    Create a SIF (simple interaction file) format of a network. Network nodes
    are reactions. SIF is compatible with Cytoscape. The filepath provided
    will be overwritten.

    :param network: Network to produce file for
    :type network: PyFBA.network.Network
    :param filepath: Filepath for output
    :type filepath: str
    :return: None
    """
    fh = open(filepath + ".sif", "w")
    conns = set()  # Record which connections were already passed through
    counts = dict()  # Record how often a reaction had a connection
    # Iterate through nodes and their neighbors
    # Algorithm: for each node, obtain all reactions connecting to it. The node
    # is essentially the connection between two reactions
    for n, nbrs in network.node_adj_iter():
        # Skip if the number of neighbors is less than two
        if len(nbrs) < 2:
            continue
        # Obtain all combinations of edges between the current
        # node and its neighboring nodes
        for n1, n2 in combinations(nbrs, 2):
            r1 = nbrs[n1]['reaction']
            r2 = nbrs[n2]['reaction']
            # Skip connections between the same reaction
            if r1 == r2:
                continue
            # Skip connections we've already seen
            check = (r1, r2)
            if check in conns:
                continue
            if r1 not in counts:
                counts[r1] = 0
            if r2 not in counts:
                counts[r2] = 0
            counts[r1] += 1
            counts[r2] += 1
            conns.add((r1, r2))
            conns.add((r2, r1))
            fh.write(r1 + "\trr\t" + r2 + "\n")
    fh.close()

    with open(filepath + ".sizes", "w") as fh:
        for rxn, count in counts.items():
            fh.write(rxn + "\t" + str(count) + "\n")


def compare_networks_sif(network1, network2, filepath, ulimit=None):
    """

    :param network1: Network1 to produce file for
    :type network1: PyFBA.network.Network
    :param network2: Network2 to produce file for
    :type network2: PyFBA.network.Network
    :param filepath: Filepath for output
    :type filepath: str
    :param ulimit: For a compound, the highest allowed number of neighbors
                   to be included in the output
    :type ulimit: int
    :return: None
    """