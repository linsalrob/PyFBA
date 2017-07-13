from __future__ import print_function, division, absolute_import
import sys
import os.path
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
        for n in self.graph.nodes():
            print(n)
        print("Number nodes before edges: " + str(len(self.graph)))

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


def clustering_coeff(network):
    pass
