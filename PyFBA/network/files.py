from __future__ import print_function, absolute_import, division
import sys
from itertools import combinations
from .metrics import network_union


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
        with open(filepath + ".sif", "w") as fh:
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

        with open(filepath + ".sizes", "w") as fh:
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

        with open(filepath + ".sif", "w") as fh:
            for c in conns_to_output:
                cpd1, cpd2 = c
                if counts[cpd1] > ulimit or counts[cpd2] > ulimit:
                    continue
                fh.write(cpd1.name + "\tcc\t" + cpd2.name + "\n")

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


def union_networks_sif(network1, network2, n1_name, n2_name,
                       filepath, ulimit=None):
    """
    Create a SIF (simple interaction file) format of the union of two networks.
    Network nodes are compounds and interactions are reactions. SIF is
    compatible with Cytoscape. The filepath provided will be overwritten.

    :param network1: Network1 to produce file for
    :type network1: PyFBA.network.Network
    :param network2: Network2 to produce file for
    :type network2: PyFBA.network.Network
    :param n1_name: Network1 name
    :type n1_name: str
    :param n2_name: Network2 name
    :type n2_name: str
    :param filepath: Filepath for output
    :type filepath: str
    :param ulimit: For a compound, the highest allowed number of neighbors
                   to be included in the output
    :type ulimit: int
    :return: None
    """
    # Obtain union of edges between two networks
    union_net = network_union(network1, network2)

    # Make networks undirected
    union_graph = union_net.get_nx_graph().to_undirected()
    g1 = network1.get_nx_graph().to_undirected()
    g2 = network2.get_nx_graph().to_undirected()

    conns = set()
    conns_to_output = set()
    counts = dict()
    node_owners = dict()

    if ulimit is None:
        pass

    else:
        # Iterate through edges and label as being unique or shared
        for e in union_graph.edges_iter():
            # Skip connection if already seen
            if e in conns:
                continue
            cpd1, cpd2 = e
            conns.add(e)
            conns.add((cpd2, cpd1))

            # Label edge as network 1, network 2, or shared
            if g1.has_edge(*e):
                if g2.has_edge(*e):
                    owner = "shared"
                else:
                    owner = n1_name
            else:
                owner = n2_name

            # Label node as network 1, network 2, or shared
            if g1.has_node(cpd1):
                if g2.has_node(cpd1):
                    node_owners[cpd1] = "shared"
                else:
                    node_owners[cpd1] = n1_name
            else:
                node_owners[cpd1] = n2_name
            if g1.has_node(cpd2):
                if g2.has_node(cpd2):
                    node_owners[cpd2] = "shared"
                else:
                    node_owners[cpd2] = n1_name
            else:
                node_owners[cpd2] = n2_name

            if cpd1 not in counts:
                counts[cpd1] = 0
            if cpd2 not in counts:
                counts[cpd2] = 0
            counts[cpd1] += 1
            counts[cpd2] += 1
            conns_to_output.add((cpd1, cpd2, owner))

        with open(filepath + ".sif", "w") as fh, \
                open(filepath + ".nodes", "w") as fh_nodes, \
                open(filepath + ".sizes", "w") as fh_sizes, \
                open(filepath + ".skip", "w") as fh_skip:
            # fh_nodes will have node assignments for each network
            fh_nodes.write("compound\tnetwork\n")
            # fh_sizes will have number of times the compound was found
            fh_sizes.write("compound\tcount\n")
            for c in conns_to_output:
                cpd1, cpd2, owner = c
                cnt1 = counts[cpd1]
                cnt2 = counts[cpd2]
                if cnt1 > ulimit:
                    fh_skip.write(cpd1.name + "\t" + str(cnt1) + "\n")
                    continue
                elif cnt2 > ulimit:
                    fh_skip.write(cpd2.name + "\t" + str(cnt2) + "\n")
                    continue
                fh.write("\t".join([cpd1.name, owner, cpd2.name]) + "\n")
                fh_nodes.write(cpd1.name + "\t" + node_owners[cpd1] + "\n")
                fh_nodes.write(cpd2.name + "\t" + node_owners[cpd2] + "\n")
                fh_sizes.write(cpd1.name + "\t" + str(cnt1) + "\n")
                fh_sizes.write(cpd2.name + "\t" + str(cnt2) + "\n")
