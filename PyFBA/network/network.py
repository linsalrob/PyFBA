from __future__ import print_function, division, absolute_import
import networkx as nx
from copy import deepcopy
import PyFBA


class Network:
    """
    A network is a directed graph of compounds (nodes) and reactions (edges).

    A network will provide methods for calculating several network metrics.

    :ivar graph: NetworkX directed graph
    """

    def __init__(self, model=None, graph=None):
        """
        Initiate the object. model takes precendence over graph

        :param model: Model object to build graph form
        :type model: PyFBA.Model
        :param graph: Graph to build Network from
        :type graph: networkx.Graph
        """
        if model is not None:
            # Initiate networkx graph
            self.graph = nx.DiGraph()
            self.__build_graph(model)
            print("Network from model " + model.name + " created!")

        elif graph is not None:
            # Check if the graph is not directed
            # If not, change to directed
            # Note: this forces all edges to form, making each bi-directional
            if not graph.is_directed():
                self.graph = deepcopy(graph.to_directed())
            else:
                self.graph = deepcopy(graph)

        else:
            raise ValueError("Model or graph is required")

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
                        # Check if edge already exists so we can append
                        # the new reaction ID
                        if self.graph.has_edge(cl, cr):
                            data = self.graph.get_edge_data(cl, cr)
                            data["reaction"] += ";" + str(rname)
                        else:
                            data = {"reaction": str(rname)}
                        self.graph.add_edge(cl, cr, data)

                    elif r.direction == "<":
                        # Check if edge already exists so we can append
                        # the new reaction ID
                        if self.graph.has_edge(cr, cl):
                            data = self.graph.get_edge_data(cr, cl)
                            data["reaction"] += ";" + str(rname)
                        else:
                            data = {"reaction": str(rname)}
                        self.graph.add_edge(cr, cl, data)

                    else:
                        # First edge
                        # Check if edge already exists so we can append
                        # the new reaction ID
                        if self.graph.has_edge(cl, cr):
                            data = self.graph.get_edge_data(cl, cr)
                            data["reaction"] += ";" + str(rname)
                        else:
                            data = {"reaction": str(rname)}

                        self.graph.add_edge(cl, cr, data)

                        # Second edge
                        if self.graph.has_edge(cr, cl):
                            data = self.graph.get_edge_data(cr, cl)
                            data["reaction"] += ";" + str(rname)
                        else:
                            data = {"reaction": str(rname)}

                        self.graph.add_edge(cr, cl, data)

    def get_nx_graph(self):
        """
        Getter for graph
        :return: Graph
        :rtype: networkx.Graph
        """
        return self.graph

    @staticmethod
    def common_compounds():
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

    def has_edge(self, edge):
        """
        Check if edge exists in network

        :param edge: Edge to check
        :type edge: tuple of PyFBA.Compound
        :return: If edge exists
        :rtype: bool
        """
        return self.graph.has_edge(*edge)

    def has_node(self, node):
        """
        Check if node exists in network

        :param node: Node to check
        :type node: PyFBA.Compound
        :return: If node exists
        :rtype: bool
        """
        return self.graph.has_node(node)

    def nodes_iter(self):
        """
        Provide an iterator for the network's nodes

        :return: Node iterator
        :rtype: iterator
        """
        for n in self.graph.nodes():
            yield n

    def nodes_all_neighbors_iter(self, node):
        """
        Provide an iterator for a node's neighbors, both predecessors and
        successors

        :param node: Node to provide an iterator for
        :type node: PyFBA.Compound
        :return: Neighbor iterator
        """
        # First make sure the node exists in the network
        if not self.graph.has_node(node):
            raise ValueError("Node does not exist in network")

        # Yield all neighboring nodes
        for n in nx.all_neighbors(self.graph, node):
            yield n

    def nodes_neighbors_iter(self, node):
        """
        Provide an iterator for a node's neighbors, only successors

        :param node: Node to provide an iterator for
        :type node: PyFBA.Compound
        :return: Neighbor iterator
        """
        # First make sure the node exists in the network
        if not self.graph.has_node(node):
            raise ValueError("Node does not exist in network")

        # Yield all successor nodes
        for n in self.graph.neighbors_iter(node):
            yield n

    def number_neighbors(self, node, all=False):
        """
        Sum number of neighbors for a node

        :param node: Node given
        :type node: PyFBA.Compound
        :param all: Flag to include predecessors in the count
        :type all: bool
        :return: Number of neighbors
        :rtype: int
        """
        # First make sure the node exists in the network
        if not self.graph.has_node(node):
            raise ValueError("Node does not exist in network")

        if all:
            return sum(1 for n in self.nodes_all_neighbors_iter(node))

        else:
            return sum(1 for n in self.nodes_neighbors_iter(node))

    def node_adj_iter(self):
        """
        Provide an iterator for each node and their adjacent nodes

        :return: Adjacency iterator
        :rtype: iterator
        """
        for n, nbrs in self.graph.adjacency_iter():
            yield (n, nbrs)

    def edges_iter(self, data=False):
        """
        Provide an iterator for the network's edges

        :param data: Flag to include edge data
        :type data: bool
        :return: Edges iterator
        :rtype: iterator
        """
        for e in self.graph.edges(data=data):
            yield e
