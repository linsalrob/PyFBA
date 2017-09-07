from __future__ import print_function, division, absolute_import
import networkx as nx
import PyFBA
import sys
from itertools import combinations

class Network:
    """
    A network is a directed graph of compounds (nodes) and reactions (edges).

    A network will provide methods for calculating several network metrics.

    :ivar graph: NetworkX directed graph
    """

    def __init__(self, model):
        """
        Initiate the object.

        :param model: Model object to build graph form
        :type model: PyFBA.Model
        """
        # Initiate networkx graph
        self.graph = nx.DiGraph()
        # self.reaction_edges = dict()

        self.__build_graph(model)
        print("Network from model " + model.name + " created!")
        print("Network contains {} nodes ".format(len(self.graph)), end="")
        print("and {} edges".format(self.graph.size()))

    def __build_graph(self, model):
        """
        Build the initial graph from the Model object.

        :param model: Model object to build graph from
        :type model: PyFBA.Model
        :return: None
        """
        # Iterate through compounds and generate nodes
        for c in model.compounds:
            self.graph.add_node(c)

        # Iterate through reactions and add edges from each reactant compound
        # to each product compound
        for r in model.reactions.values():
            # Instantiate reaction edge in dict
            rname = str(r)
            # self.reaction_edges[rname] = [set(), set(), r.direction]
            for cl in r.left_compounds:
                # Add left compound to edge dict
                # self.reaction_edges[rname][0].add(str(cl))
                if not self.graph.has_node(cl):
                    print("New compound: {}, reaction: {}".format(cl, r))
                for cr in r.right_compounds:
                    # Add right compound to edge dict
                    # self.reaction_edges[rname][1].add(str(cr))
                    if not self.graph.has_node(cr):
                        print("New compound: {}, reaction: {}".format(cr, r))
                    # Since graph is directed, the order of the compounds
                    # in the add_edge() function matters
                    if r.direction == ">":
                        self.graph.add_edge(cl, cr, {"reaction" : str(rname)})

                    elif r.direction == "<":
                        self.graph.add_edge(cr, cl, {"reaction" : str(rname)})

                    else:
                        self.graph.add_edge(cl, cr, {"reaction" : str(rname)})
                        self.graph.add_edge(cr, cl, {"reaction" : str(rname)})

    def common_compounds(self):
        """
        Return set of highly common compounds
        """
        cc = ["H+", "H2O", "ATP", "ADP", "Phosphate",
              "NAD", "NADH", "NADP", "NADPH", "PPi",
              "CoA", "CO2"]
        cpds = set()
        for c in cc:
            cpds.add(PyFBA.metabolism.Compound(c, "c"))
            cpds.add(PyFBA.metabolism.Compound(c, "e"))
        return cpds

    def number_of_nodes(self):
        """
        Provide number of nodes (compounds) in the network

        :return: Number of nodes
        :rtype: int
        """
        return self.graph.number_of_nodes()

    def number_of_edges(self):
        """
        Provide number of edges (reactions) in the network

        :return: Number of edges
        :rtype: int
        """
        return self.graph.number_of_edges()

    def nodes_iter(self):
        """
        Provide an iterator for the network's nodes

        :return: Node iterator
        :rtype: iterator
        """
        for n in self.graph.nodes():
            yield n

    def nodes_neighbors_iter(self, node):
        """
        Provide an iterator for a node's neighbors

        :param node: Node to provide an iterator for
        :type node: PyFBA.Compound
        :return: Neighbor iterator
        """
        # First make sure the node exists in the network
        if not self.graph.has_node(node):
            raise ValueError("Node does not exist in network")

        # Make network undirected first
        for n in self.graph.to_undirected().neighbors_iter(node):
            yield n

    def number_neighbors(self, node):
        """
        Sum number of neighbors for a node

        :param node: Node given
        :type node: PyFBA.Compound
        :return: Number of neighbors
        :rtype: int
        """
        # First make sure the node exists in the network
        if not self.graph.has_node(node):
            raise ValueError("Node does not exist in network")

        # Make network undirected first
        return sum(1 for n in self.nodes_neighbors_iter(node))

    def node_adj_iter(self):
        """
        Provide an iterator for each node and their adjacent nodes

        :return: Adjacency iterator
        :rtype: iterator
        """
        for n, nbrs in self.graph.to_undirected().adjacency_iter():
            yield (n, nbrs)

    def edges_iter(self):
        """
        Provide an iterator for the network's edges

        :return: Edges iterator
        :rtype: iterator
        """
        for e in self.graph.edges():
            yield e


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
        counts = dict()
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
        for cpd, count in counts.items():
            fhdata.write(cpd.name + "\t" + count + "\n")
        fh.close()
        fhdata.close()

    else:
        counts = dict()
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
