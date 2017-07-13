from __future__ import print_function, absolute_import, division
import sys
import os
import PyFBA
from .network import Network
import networkx as nx


def clustering_coeff(network):
    """
    Calculate the clustering coefficient using the NetworkX library.

    :param network: Network of nodes
    :type network: PyFBA.network.Network
    :return: Clustering coefficient for nodes in the network
    :rtype: dict
    """
    # Check type
    if not isinstance(network, Network):
        raise TypeError

    # Check if network is empty
    if len(network.graph) == 0:
        raise ValueError("Network is empty")

    # Convert directed graph to undirected
    if network.graph.is_directed():
        return nx.clustering(network.graph.to_undirected())

    return nx.clustering(network.graph)


def avg_clustering_coeff(network):
    """
    Calculate the average clustering coefficient using the NetworkX library.

    :param network: Network of nodes and edges
    :type network: PyFBA.network.Network
    :return: Average clustering coefficient for all nodes in the network
    :rtype: float
    """
    # Check type
    if not isinstance(network, Network):
        raise TypeError

    # Check if network is empty
    if len(network.graph) == 0:
        raise ValueError("Network is empty")

    return nx.average_clustering(network.graph)


def closeness_centrality(network):
    """
    Calculate the closeness centrality using the NetworkX library.

    :param network: Network of nodes and edges
    :type network: PyFBA.network.Network
    :return: Closeness centrality for nodes in the network
    :rtype: dict
    """
    # Check type
    if not isinstance(network, Network):
        raise TypeError

    # Check if network is empty
    if len(network.graph) == 0:
        raise ValueError("Network is empty")

    return nx.closeness_centrality(network.graph)


def avg_closeness_centrality(network):
    """
    Calculate the average closeness centrality using the NetworkX library.

    :param network: Network of nodes and edges
    :type network: PyFBA.network.Network
    :return: Average closeness centrality for all nodes in the network
    :rtype: float
    """
    cc = closeness_centrality(network)
    return sum(list(cc.values())) / len(cc)


def betweenness_centrality(network):
    """
    Calculate the betweenness centrality using the NetworkX library.

    :param network: Network of nodes and edges
    :type network: PyFBA.network.Network
    :return: Betweenness centrality for nodes in the network
    :rtype: dict
    """
    # Check type
    if not isinstance(network, Network):
        raise TypeError

    # Check if network is empty
    if len(network.graph) == 0:
        raise ValueError("Network is empty")

    return nx.betweenness_centrality(network.graph)


def avg_betweenness_centrality(network):
    """
    Calculate the average betweenness centrality using the NetworkX library.

    :param network: Network of nodes and edges
    :type network: PyFBA.network.Network
    :return: Average betweenness centrality for all nodes in the network
    :rtype: float
    """
    bc = betweenness_centrality(network)
    return sum(list(bc.values())) / len(bc)


def diameter(network):
    """
    Calculate the diameter of the network (the longest path between two nodes).

    :param network: Network of nodes and edges
    :type network: PyFBA.network.Network
    :return: Network diameter
    :rtype: int
    """
    # Check type
    if not isinstance(network, Network):
        raise TypeError

    # Check if network is empty
    if len(network.graph) == 0:
        raise ValueError("Network is empty")

    return nx.diameter(network.graph)
