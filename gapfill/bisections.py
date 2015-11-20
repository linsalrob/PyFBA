def bisect(rxns):
    """
    Given a list of reactions, split it in two. This just returns two
    lists, with 1/2 the elements each. Note that this is currently an in
    order split but that the lists are elements [0,2,4,6...] and
    [1,3,5,7...]

    :param rxns: A list of reactions
    :type rxns: list
    :return: Two lists, evenly split
    :rtype: list, list
    """

    rlist = list(rxns)
    # [::2] gets every 2nd element. 1::2 does it from the 2nd element
    return rlist[::2], rlist[1::2]


def percent_split(rxns, percent=50):
    """
    Given a list of reactions, split the list into two different lists,
    one of which will be a percentage of the total, and the other
    will be the remainder. 

    For example, if percentage is 20, one list will have 20% of the
    elements and the other list will have 80% of the elements.

    IF the percentage is 50 (default), each list has 50% of the elements

    :param percent: The percentage in which to split the list
    :type percent: int
    :param rxns: A list of reactions
    :type rxns: list
    :return: Two lists, evenly split
    :rtype: list, list
    """

    brk = int(len(rxns) * (1.0 * percent/100))
    return rxns[:brk], rxns[brk:]


def optimize_split_by_rclust(rxns, clusters, percent=50):
    """
    Optimize the split of reactions based on the reaction clusters.

    We have code that clusters reactions based on the Jaccard distance between the compounds in the reactions (see
    example_code/group_reactions.py). This code uses that cluster information to try and optimize the split in reactions
    to minimize the number of different clusters on each side of the split.

    The approach that we take here is to sort the reactions based on their cluster first, and then split the list
    although this imposes an order on the list before sort it should still work (I think) while being a lot more
    efficient than sorting after the split (see post_optimize_split_by_rclust)

    :param percent: The percent of the split: 50% is an even split
    :type percent: int
    :param rxns: The list of reactions to split
    :type rxns: list
    :param clusters: a hash with the reaction ids and the clusters they are in
    :type clusters: dict
    :return: the two lists of reactions
    :rtype: (list, list)
    """

    rxnssort = sorted(rxns, key=clusters.get)
    return percent_split(rxnssort, percent)

