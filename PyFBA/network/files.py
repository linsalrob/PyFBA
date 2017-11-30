from __future__ import print_function, absolute_import, division
import sys
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
    counts = dict()

    num_nodes = network.number_of_nodes()
    num_edges = network.number_of_edges()
    if verbose:
        print("Network contains", num_nodes, "nodes and", num_edges, "edges",
              file=sys.stderr)

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
        if cpd1 not in counts:
            counts[cpd1] = 0
        if cpd2 not in counts:
            counts[cpd2] = 0
        counts[cpd1] += 1
        counts[cpd2] += 1
    if verbose:
        print("", file=sys.stderr)

    # Write out files
    with open(filepath + ".sif", "w") as fh, \
            open(filepath + ".sizes", "w") as fh_sizes, \
            open(filepath + ".skip", "w") as fh_skip:
        # Remember which compounds are seen to avoid repetition in output files
        reported_skip = set()
        reported_sizes = set()
        for c in conns:
            skip = False
            cpd1, cpd2 = c
            cnt1 = counts[cpd1]
            cnt2 = counts[cpd2]
            cpd1_name = "{} [{}]".format(cpd1.name, cpd1.location)
            cpd2_name = "{} [{}]".format(cpd2.name, cpd2.location)

            # Check if first compound is too high
            if ulimit and cnt1 > ulimit:
                # Check if first compound was written to file yet
                if cpd1 not in reported_skip:
                    fh_skip.write(cpd1_name + "\t" + str(cnt1) + "\n")
                    reported_skip.add(cpd1)
                skip = True

            # Check if second compound is too high
            if ulimit and cnt2 > ulimit:
                # Check if second compound was written to file yet
                if cpd2 not in reported_skip:
                    fh_skip.write(cpd2_name + "\t" + str(cnt2) + "\n")
                    reported_skip.add(cpd2)
                skip = True
            if skip:
                continue

            # Check if compounds were written to sizes file yet
            if cpd1 not in reported_sizes:
                fh_sizes.write(cpd1_name + "\t" + str(cnt1) + "\n")
                reported_sizes.add(cpd1)
            if cpd2 not in reported_sizes:
                fh_sizes.write(cpd2_name + "\t" + str(cnt2) + "\n")
                reported_sizes.add(cpd2)

            # Write out to interaction file
            fh.write(cpd1_name + "\tcc\t" + cpd2_name + "\n")


def network_reactions_sif(network, filepath, ulimit=None):
    """
    Create a SIF (simple interaction file) format of a network. Network nodes
    are reactions. SIF is compatible with Cytoscape. The filepath provided
    will be overwritten.

    :param network: Network to produce file for
    :type network: PyFBA.network.Network
    :param filepath: Filepath for output
    :type filepath: str
    :param ulimit: For a compound, the highest allowed number of neighbors
                   to be included in the output
    :type ulimit: int
    :return: None
    """
    # Record for each compound the reactions in which the
    # compound is a reactant and the reactions in which the
    # compound is a product
    cpds = dict()

    conns = set()  # Record which connections were already passed through
    counts = dict()  # Record how often a reaction had a connection

    # Iterate through nodes (compounds) and their neighbors
    # Algorithm: for each node (compound), obtain all reactions connecting to
    # it. The node is essentially the connection between two reactions
    for c, nbrs in network.node_adj_iter():
        # Neighbors stored as:
        #     {cpd_neighbor1: {"reaction": rxn_id},
        #      cpd_neighbor2: {"reaction": rxn_id}, ...}
        for n, data in nbrs.items():
            r = data["reaction"]
            # Add right reaction to the current compound
            if c not in cpds:
                cpds[c] = {"left": set(), "right": set()}
            cpds[c]["right"].add(r)

            # Add the left reaction to the neighboring compounds
            if n not in cpds:
                cpds[n] = {"left": set(), "right": set()}
            cpds[n]["left"].add(r)

    # Iterate through all compounds
    # Reactions are connected by the left of the compound and the right
    # of the compound
    for c in cpds:
        for r1 in cpds[c]["left"]:
            for r2 in cpds[c]["right"]:
                # Skip connections we've already seen
                # This shouldn't happen anymore but we'll keep it here anyway
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

    # Write out files
    with open(filepath + ".sif", "w") as fh, \
            open(filepath + ".sizes", "w") as fh_sizes, \
            open(filepath + ".skip", "w") as fh_skip:
        # Remember which compounds are seen to avoid repetition in output files
        reported_skip = set()
        reported_sizes = set()
        for c in conns:
            skip = False
            rxn1, rxn2 = c
            cnt1 = counts[rxn1]
            cnt2 = counts[rxn2]

            # Check if first compound is too high
            if ulimit and cnt1 > ulimit:
                # Check if first compound was written to file yet
                if rxn1 not in reported_skip:
                    fh_skip.write(rxn1 + "\t" + str(cnt1) + "\n")
                    reported_skip.add(rxn1)
                skip = True

            # Check if second compound is too high
            if ulimit and cnt2 > ulimit:
                # Check if second compound was written to file yet
                if rxn2 not in reported_skip:
                    fh_skip.write(rxn2 + "\t" + str(cnt2) + "\n")
                    reported_skip.add(rxn2)
                skip = True
            if skip:
                continue

            # Check if compounds were written to sizes file yet
            if rxn1 not in reported_sizes:
                fh_sizes.write(rxn1 + "\t" + str(cnt1) + "\n")
                reported_sizes.add(rxn1)
            if rxn2 not in reported_sizes:
                fh_sizes.write(rxn2 + "\t" + str(cnt2) + "\n")
                reported_sizes.add(rxn2)

            # Write out to interaction file
            fh.write(rxn1 + "\trr\t" + rxn2 + "\n")


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
    union_graph = union_net.get_nx_graph()

    conns = set()
    conns_to_output = set()
    counts = dict()
    node_owners = dict()

    # Iterate through edges and label as being unique or shared
    for e in union_graph.edges_iter():
        # Skip connection if already seen
        if e in conns:
            continue
        cpd1, cpd2 = e
        conns.add(e)
        conns.add((cpd2, cpd1))

        # Label edge as network 1, network 2, or shared
        if network1.has_edge((cpd1, cpd2)):
            if network2.has_edge((cpd1, cpd2)):
                owner = "shared"
            else:
                owner = n1_name
        else:
            owner = n2_name

        # Label node as network 1, network 2, or shared
        if network1.has_node(cpd1):
            if network2.has_node(cpd1):
                node_owners[cpd1] = "shared"
            else:
                node_owners[cpd1] = n1_name
        else:
            node_owners[cpd1] = n2_name
        if network1.has_node(cpd2):
            if network2.has_node(cpd2):
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
        # Remember which compounds are seen to avoid repetition in
        # output files
        reported_skip = set()
        reported_others = set()
        for c in conns_to_output:
            skip = False
            cpd1, cpd2, owner = c

            cnt1 = counts[cpd1]
            cnt2 = counts[cpd2]
            cpd1_name = "{} [{}]".format(cpd1.name, cpd1.location)
            cpd2_name = "{} [{}]".format(cpd2.name, cpd2.location)

            # Check if first compound count is too high
            if ulimit is not None and cnt1 > ulimit:
                # Check if first compound was written to file yet
                if cpd1 not in reported_skip:
                    fh_skip.write(cpd1_name + "\t" + str(cnt1) + "\n")
                    reported_skip.add(cpd1)
                skip = True

            # Check if second compound count is too high
            if ulimit is not None and cnt2 > ulimit:
                # Check if second compound was written to file yet
                if cpd2 not in reported_skip:
                    fh_skip.write(cpd2_name + "\t" + str(cnt2) + "\n")
                    reported_skip.add(cpd2)
                skip = True
            if skip:
                continue

            # Check if compounds were written to nodes and sizes file yet
            if cpd1 not in reported_others:
                fh_nodes.write(cpd1_name + "\t" + node_owners[cpd1] + "\n")
                fh_sizes.write(cpd1_name + "\t" + str(cnt1) + "\n")
                reported_others.add(cpd1)
            if cpd2 not in reported_others:
                fh_nodes.write(cpd2_name + "\t" + node_owners[cpd2] + "\n")
                fh_sizes.write(cpd2_name + "\t" + str(cnt2) + "\n")
                reported_others.add(cpd2)

            # Write out to interaction file
            fh.write("\t".join([cpd1_name, owner, cpd2_name]) + "\n")
