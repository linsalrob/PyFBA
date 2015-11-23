"""
Given a set of reactions (e.g. all of our biochemistry) can we come up with some sane groups of reactions?
This should segregate reactions based on common compounds, but not those that are present in everything (e.g. H2O)
"""
import argparse
import sys

import PyFBA


def jaccard(s1, s2):
    """
    Calculate the Jaccard *distance* between two sets (the length of the intersection / the length of the union). If
    either or both sets are empty the distance is 1.

    Note that this returns a distance -- 0 means things are similar, 1 means things are different

    :param s1: Set one
    :type s1: set
    :param s2: Set two
    :type s2: set
    :return: The Jaccard distance
    :rtype: float
    """

    if len(s1) == 0 or len(s2) == 0:
        return 1
    return 1 - 1.0 * len(s1.intersection(s2)) / len(s1.union(s2))


def calculate_distances(reactions, threshold):
    """
    Calculate the pairwise Jaccard distances for all pairs of reactions
    :param reactions: the reactions dictionary
    :type reactions: dict
    :param threshold: The threshold for which to join clusters
    :type threshold: float
    :return: A hash of a reaction, another reaction, and the distance between them
    :rtype: dict
    """
    # Calculate a pairwise distance between all the reactions in our biochemistry
    distance = {}
    rcts = reactions.keys()
    for s in rcts:
        if s not in distance:
            distance[s] = {}
        for t in rcts:
            dist = jaccard(reactions[s].all_compounds(), reactions[t].all_compounds())
            if dist > threshold:
                continue
            if t in distance and s in distance[t]:
                continue
            if t not in distance:
                distance[t] = {}
            # the distance is the number of differences they have
            distance[s][t] = distance[t][s] = dist

    return distance


def write_distances(reactions, filename):
    """
    Write the distances to a file (filename) rather than storing in memory. You may need to use this if you
    don't have a lot of RAM
    :param reactions: The reactions dictionary
    :type reactions: dict
    :param filename: The filename to use to store the distances
    :type filename: str
    :return: Nothing
    :rtype:
    """
    # Calculate a pairwise distance between all the reactions in our biochemistry
    rcts = reactions.keys()

    with open(filename, 'w') as out:
        for s in rcts:
            for t in rcts:
                # the distance is the number of differences they have
                out.write("\t".join([s, t, str(jaccard(reactions[s].all_compounds(),
                                                       reactions[t].all_compounds()))]) + "\n")


def calculate_clusters(reactions, distance, threshold):
    """
    Calculate the clusters of reactions based on their distances
    :param reactions: A list of all the reactions
    :type reactions: list
    :param distance: The hash of reaction/reaction and distances
    :type distance: dict
    :param threshold: The threshold to consider clustering
    :type threshold: float
    :return: A hash of the reactions and the clusters to which they belong and the size of the clusters
    :rtype: (dict, int)
    """
    # Cluster the distances using single linkage clustering
    cluster = {}
    current_cluster = 0
    for r in reactions:
        # first test if this reaction is linked to an existing cluster
        if r not in distance:
            cluster[r] = current_cluster
            current_cluster += 1
            continue
        for s in cluster:
            if s not in distance[r]:
                continue
            if distance[r][s] < threshold:
                cluster[r] = cluster[s]
                break
        if r not in cluster:
            cluster[r] = current_cluster
            current_cluster += 1

    return cluster, current_cluster - 1


def read_distance_calculate_clusters(filename, reactions, threshold):
    """
    Calculate the clusters of reactions based on their distances stored in a file
    :param reactions: A list of all the reactions in the model
    :type reactions: list
    :param filename: The filename that has [rxn1, rxn2, distance]
    :type filename: str
    :param threshold: The threshold to use for joining a cluster
    :type threshold: float
    :return: A hash of the reaction ids and their clusters and the number of clusters
    :rtype: (dict, int)
    """

    cluster = {}
    rxnbycluster = {}
    current_cluster = 0
    with open(filename, 'r') as f:
        for l in f:
            rfrom, rto, dist = l.strip().split("\t")
            if dist > threshold:
                continue
            if rfrom in cluster and rto in cluster and cluster[rfrom] == cluster[rto]:
                # in the same cluster, good
                continue
            elif rfrom in cluster and rto in cluster:
                # these clusters should be merged
                minclust = min(cluster[rfrom], cluster[rto])
                maxclust = max(cluster[rfrom], cluster[rto])
                for r in rxnbycluster[maxclust]:
                    cluster[r] = minclust
                # remove maxclust from rxnbycluster and also update the set for min
                # this should ensure we get an error if we try and access maxclust again!
                rxnbycluster[minclust].update(rxnbycluster.pop(maxclust))
                continue
            elif rfrom in cluster:
                cluster[rto] = cluster[rfrom]
                rxnbycluster[cluster[rfrom]].add(rto)
            elif rto in cluster:
                cluster[rfrom] = cluster[rto]
                rxnbycluster[cluster[rto]].add(rfrom)
            else:
                cluster[rfrom] = current_cluster
                cluster[rto] = current_cluster
                rxnbycluster[current_cluster] = {rfrom, rto}
                current_cluster += 1

    # now add the missing reactions
    for r in reactions:
        if r not in cluster:
            cluster[r] = current_cluster
            current_cluster += 1

    return cluster, current_cluster - 1


def group_reactions(reactions, rcts, compound_threshold, verbose=False):
    """
    Group the reactions based on the connectivity of compounds. Only compounds with a distance < threshold are
    considered and we use greedy clustering (ie. if any compound is in a group we add it).

    The distance is calculated using the Jaccard distance (appropriate for two sets) - the length of the intersection/
    the length of the union

    :param rcts: A list of the reactions in reactions
    :type rcts: list
    :param reactions: the reactions dictionary
    :type reactions: dict
    :param compound_threshold: The threshold for the number of reactions a compound is in
    :type compound_threshold: int
    :param verbose: Whether to print out more stuff
    :type verbose: bool
    :return: A hash of rxn ids and the group they are in and the number of clusters
    :rtype: (dict, int)
    """
    d = calculate_distances(reactions, verbose)
    return calculate_clusters(rcts, d, compound_threshold)


def read_distance_file(dist_file, threshold):
    """
    Read a previously created distance file and store it as a hash
    :param threshold: The threshold for the distances
    :type threshold: float
    :param dist_file: The file to read
    :type dist_file: str
    :return: A hash of rxn1 rxn2 -> distance
    :rtype: dict
    """
    distances = {}
    with open(dist_file, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if float(p[2]) > threshold:
                continue
            if p[0] not in distances:
                distances[p[0]] = {}
            if p[1] not in distances:
                distances[p[1]] = {}
            distances[p[0]][p[1]] = float(p[2])
            distances[p[1]][p[0]] = float(p[2])
    return distances


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculate and write out clusters of reactions based on compound similarity')
    parser.add_argument('-f', help="Use file as interim file to store distances (uses a lot less memory!")
    parser.add_argument('-d',
                        help="Distances file. If you have this computed you can iterate through the clusters to " +
                             "explore the sizes (required with -c)")
    parser.add_argument('-c', help="Count the number of reactions in clusters with different threshold sizes",
                        action='store_true')
    parser.add_argument('-t', help="Threshold for clustering. Default = 0.5", type=float, default=0.5)
    args = parser.parse_args()

    compounds, allreactions, enzymes = PyFBA.parse.model_seed.compounds_reactions_enzymes()

    reaction_keys = allreactions.keys()

    if args.f:
        write_distances(allreactions, args.f)
        clusters = read_distance_calculate_clusters(args.f, reaction_keys, args.t)
        for c in clusters:
            print("\t".join([c, str(clusters[c])]))
    elif args.c:
        if not args.d:
            sys.exit(
                "The -d parameter must be defined to work with the -c option. Please run first and create that file!")
        sys.stderr.write("Threshold | Number of clusters\n --- | ---\n")
        for i in range(1, 10):
            th = 1.0 * i / 10
            dists = read_distance_file(args.d, th)
            clusters, cluster_len = calculate_clusters(reaction_keys, dists, th)
            sys.stderr.write(str(th) + " | " + str(cluster_len) + "\n")
            for c in clusters:
                print('\t'.join([str(th), c, str(clusters[c])]))
    elif args.d:
        dists = read_distance_file(args.d, args.t)
        clusters, cluster_len = calculate_clusters(reaction_keys, dists, args.t)
        sys.stderr.write(
            "File: " + args.d + " threshold: " + str(args.t) + " number of clusters: " + str(cluster_len) + "\n")
        for c in clusters:
            print("\t".join([args.d, c, str(clusters[c])]))
    else:
        clusters = group_reactions(allreactions, reaction_keys, args.t)
        for c in clusters:
            print("\t".join([c, str(clusters[c])]))
