import copy
import sys
from random import shuffle

import time
from gapfill import bisections, limit_reactions_by_compound
import Queue
from threading import Thread

__author__ = 'Rob Edwards'


class FBAWorker(Thread):
    import fba

    def __init__(self, queue):
        Thread.__init__(self)
        self.queue = queue

    def run(self):
        while True:
            params = self.queue.get()
            status, value, growth = self.__class__.fba.run_fba(*params)
            sys.stderr.write("IN WORKER HAD {} and {}\n".format(value, growth))
            self.queue.task_done()





def minimize_additional_reactions_pl(base_reactions, optional_reactions, compounds, reactions, media,
                                  biomass_eqn, verbose=False):
    """
    Given two sets, one of base reactions (base_reactions), and one of optional
    reactions we will attempt to minimize the reactions in the optional
    reactions list by repeatedly bisecting the optional set and testing
    the fba. We return a set of the optional reactions that are required
    for the fba to grow.

    :param compounds: The compounds dictionary
    :type compounds: dict
    :param base_reactions: a set of reactions that are required for the model but that do not result in growth
    :type base_reactions: set
    :param optional_reactions: a set of reactions that when added to the base_reactions set result in
        growth but for which only a subset may or may not be required.
    :type optional_reactions: set
    :param reactions: the reactions data dictionary
    :type reactions: dict
    :param media: our media object
    :type media: set
    :param biomass_eqn: our biomass equation
    :type biomass_eqn: network.reaction.Reaction
    :param verbose: Print more information
    :type verbose: bool
    :return: The set of reactions that need to be added to base_reactions to get growth
    :rtype: set


    """

    q = Queue.Queue()

    base_reactions = set(base_reactions)
    optional_reactions = set(optional_reactions)

    # test that (a) the base_reactions set does not grow and the base_reactions + optional set does grow
    start = time.time()
    params = [
        [compounds, reactions, base_reactions, media, biomass_eqn],
        [compounds, reactions, base_reactions.union(optional_reactions), media, biomass_eqn]
    ]

    sys.stderr.write("BEFORE PASS TO FBA: c: {} r: {} r2r: {} (other r2r: {})\n".format(len(compounds),
        len(reactions), len(base_reactions), len(base_reactions.union(optional_reactions))))

    # results is a list of the results from run_fba. We only care about growth.

    #results = pool.map(pass_to_fba, params)
    for p in params:
        fbaw = FBAWorker(q)
        fbaw.daemon = True
        fbaw.start()

    for p in params:
        q.put(p)

    q.join()
    results = q.get()

    lgrowth = results[0][2]
    rgrowth = results[1][2]

    #results = fba.run_fba(*params[0])
    #lgrowth = results[2]
    #results = fba.run_fba(*params[1])
    #rgrowth = results[2]

    sys.stderr.write("Running WITH PARALLEL took {}\n".format(time.time() - start))

    if lgrowth:
        sys.stderr.write("The set of 'base' reactions results in growth so we don't need to bisect the optional set\n")
        return set()
    if not rgrowth:
        raise Exception("'base' union 'optional' reactions does not generate growth. We can not bisect the set\n")

    # now lets see if we can limit the reactions based on compounds present and still get growth
    limited_rxn = limit_reactions_by_compound(reactions, base_reactions, optional_reactions)
    # status, value, growth = fba.run_fba(compounds, reactions, base_reactions.union(limited_rxn), media,
    #                                    biomass_eqn)
    #  ############### WE NEED SOMETHING HERE!!!!!!
    # results = pass_to_fba(q, [compounds, reactions, base_reactions.union(limited_rxn), media, biomass_eqn])
    growth = results[2]
    if growth:
        if verbose:
            sys.stderr.write("Successfully limited the reactions by compound and reduced " +
                             " from {} to {}\n".format(len(optional_reactions), len(limited_rxn)))
        optional_reactions = limited_rxn

    # now iterate through our list of reactions trying to break them down to the minimal set
    test = True
    tries = 0
    maxtries = 5
    current_rx_list = list(optional_reactions)
    itera = 0
    sys.stderr.write("At the beginning the base list has {} ".format(len(base_reactions)) +
                     " and the optional list has {} reactions\n".format(len(current_rx_list)))

    left = []
    right = []
    while test:
        itera += 1
        left, right = bisections.bisect(current_rx_list)
        start = time.time()
        params = [
            [compounds, reactions, base_reactions.union(set(left)), media, biomass_eqn],
            [compounds, reactions, base_reactions.union(set(left)), media, biomass_eqn]
        ]

        # results is a list of the results from run_fba. We only care about growth.
        ############### WE NEED SOMETHING HERE!!!!!!
        lgrowth = results[0][2]
        rgrowth = results[1][2]
        sys.stderr.write("Running WITH PARALLEL took {}\n".format(time.time() - start))

        if verbose:
            sys.stderr.write("Iteration: {} Try: {} Length: {} and {}".format(itera, tries, len(left), len(right)) +
                             " Growth: {} and {}\n".format(lgrowth, rgrowth))

        if lgrowth and rgrowth:
            # they both grow, so we can choose one!
            tries = 0
            current_rx_list = left
        elif lgrowth:
            # the left list grows
            tries = 0
            current_rx_list = left
            if len(left) == 1:
                test = False
                right = []
        elif rgrowth:
            # the right list grows
            tries = 0
            current_rx_list = right
            if len(right) == 1:
                test = False
                left = []
        else:
            # neither grows. Can we split the list unevenly and see
            # if we get growth
            percent = 40
            left, right = bisections.percent_split(current_rx_list, percent)
            uneven_test = True
            while uneven_test and len(left) > 0 and len(right) > 0:
                start = time.time()
                params = [
                    [compounds, reactions, base_reactions.union(set(left)), media, biomass_eqn],
                    [compounds, reactions, base_reactions.union(set(left)), media, biomass_eqn]
                ]
                # results is a list of the results from run_fba. We only care about growth.

                ############### WE NEED SOMETHING HERE!!!!!!

                lgrowth = results[0][2]
                rgrowth = results[1][2]
                sys.stderr.write("Running WITH PARALLEL took {}\n".format(time.time() - start))
                if verbose:
                    sys.stderr.write(
                        "Iteration: {} Try: {} Length: {} and {}".format(itera, tries, len(left), len(right)) +
                        " Growth: {} and {}\n".format(lgrowth, rgrowth))
                if lgrowth:
                    tries = 0
                    current_rx_list = left
                    uneven_test = False
                elif rgrowth:
                    tries = 0
                    current_rx_list = right
                    uneven_test = False
                else:
                    percent /= 2.0
                    left, right = bisections.percent_split(current_rx_list, percent)
            if uneven_test:
                # we never got growth, so we can't continue
                # we take another stab and try again
                tries += 1
                current_rx_list = left + right
                shuffle(current_rx_list)
            if tries > maxtries:
                test = False

    remaining = set(left + right)
    if verbose:
        sys.stderr.write("There are {} reactions remaining: {}\n".format(len(remaining), remaining))
    return remaining




