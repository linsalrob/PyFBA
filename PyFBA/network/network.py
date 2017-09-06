from __future__ import print_function, division, absolute_import
import networkx as nx
import PyFBA


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
            for cl in r.left_compounds:
                if not self.graph.has_node(cl):
                    print("New compound: {}, reaction: {}".format(cl, r))
                for cr in r.right_compounds:
                    if not self.graph.has_node(cr):
                        print("New compound: {}, reaction: {}".format(cr, r))
                    # Since graph is directed, the order of the compounds
                    # in the add_edge() function matters
                    if r.direction == ">":
                        self.graph.add_edge(cl, cr)

                    elif r.direction == "<":
                        self.graph.add_edge(cr, cl)

                    else:
                        self.graph.add_edge(cl, cr)
                        self.graph.add_edge(cr, cl)

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

    def edges_iter(self):
        """
        Provide an iterator for the network's edges

        :return: Edges iterator
        :rtype: iterator
        """
        for e in self.graph.edges():
            yield e


def network_compounds_sif(network,  filepath):
    """
    Create a SIF (simple interaction file) format of a network. Network nodes
    are compounds. SIF is compatible with Cytoscape. The filepath provided
    will be overwritten.

    :param network: Network to produce file for
    :type network: PyFBA.network.Network
    :param filepath: Filepath for output
    :type filepath: str
    :return: None
    """
    fh = open(filepath, "w")
    conns = set()  # Record which connections were already passed through
    for e in network.edges_iter():
        # Skip connections we've already seen
        # Could be due to bidirectional reactions
        if e in conns:
            continue
        cpd1, cpd2 = e
        conns.add(e)
        conns.add((cpd2, cpd1))
        fh.write(cpd1.name + "\tcc\t" + cpd2.name + "\n")
    fh.close()


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
    pass
