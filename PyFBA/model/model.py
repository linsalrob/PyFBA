from __future__ import print_function
import PyFBA
import copy
import sys
import os
from os.path import basename


class Model:
    """
    A model is a structured collection of reactions, compounds, and functional roles.

    A model will provide methods for extracting its information and creating SBML-formatted files.

    :ivar id: Model ID
    :ivar name: Model name
    :ivar roles: A dictionary of roles as the key and sets of reaction IDs as the value
    :ivar reactions: A dictionary of reaction IDs as the key and Reaction objects as the value
    :ivar compounds: A set of Compound objects
    :ivar gapfilled_media: A set of media names this model has been gap-filled with
    :ivar gf_reactions: A dictionary of gap-filled reaction IDs as the key and their gap-filled step as the value
    :ivar biomass_reaction: A reaction object representing the biomass reaction
    :ivar organism_type: A String describing the type of organism
    """

    def __init__(self, id, name, organism_type="gramnegative"):
        """
        Initiate the object

        :param id: Model ID
        :type id: str
        :param name: Model name
        :type name: str
        :return:
        :rtype:
        """
        self.id = str(id)
        self.name = str(name)
        self.roles = dict()
        self.reactions = dict()
        self.compounds = set()
        self.gapfilled_media = set()
        self.gf_reactions = dict()
        self.biomass_reaction = None
        self.organism_type = organism_type

    def __str__(self):
        """
        The to string function.

        :rtype: str
        """
        return self.name + " (id: " + self.id  + ")"

    def add_roles(self, roles):
        """
        Add roles to the model.

        :param roles: A dictionary of roles and reaction IDs
        :type roles: dict of set of str
        """
        if isinstance(roles, dict):
            for role, rxns in roles.items():
                if role not in self.roles:
                    self.roles[role] = set()
                self.roles[role].update(rxns)

    def add_reactions(self, rxns):
        """
        Add reactions to the model.

        :param rxns: Reaction objects
        :type rxns: set
        """
        if isinstance(rxns, set):
            for r in rxns:
                if r.name not in self.reactions:
                    self.reactions[r.name] = r
                    # Add compounds from new reactions
                    for c in r.all_compounds():
                        if c not in self.compounds:
                            self.compounds.add(c)
        else:
            raise TypeError("You need to add a set of reactions to a model")

    def remove_reactions(self, rxns):
        """
        Remove reactions from the model.

        :param rxns: Reaction IDs
        :type rxns: set
        """
        # Passing for now
        # Need to figure out best way to deal with compounds
        pass

    def has_reaction(self, rxn):
        """
        Check if model contains reaction.

        :param rxn: A Reaction object
        :type rxn: Reaction
        :rtype: bool
        """
        return rxn.name in self.reactions

    def number_of_reactions(self):
        """
        Return number of reactions this model contains.

        :rtype: int
        """
        return len(self.reactions)

    def has_compound(self, cpd):
        """
        Check if model contains compound.

        :param cpd: A Compound object
        :type cpd: Compound
        :rtype: bool
        """
        return cpd in self.compounds

    def number_of_compounds(self):
        """
        Return number of compounds this model contains.

        :rtype: int
        """
        return len(self.compounds)

    def set_biomass_reaction(self, rxn):
        """
        Specify biomass reaction for model.

        :param rxn: A Reaction object.
        :type rxn: Reaction
        """
        self.biomass_reaction = rxn

    def output_model(self, f):
        """
        Output model reaction, function, and gap-fill information.

        :param f: File object to print to
        :type f: file
        """
        # Get mapping from reaction IDs to roles
        mReactions = {r: [] for r in self.reactions.keys()}
        for role, rxns in self.roles.items():
            for r in rxns:
                mReactions[r].append(role)
        # Print header
        f.write("reaction\tfunction\tequation\tgapfilled\n")

        # Print reactions from model
        for r, roles in mReactions.items():
            eqn = self.reactions[r].equation
            rolecolumn = ";".join(roles)
            if r in self.gf_reactions:
                f.write("\t".join([r, rolecolumn, eqn, "yes"]))
            else:
                f.write("\t".join([r, rolecolumn, eqn, "no"]))
            f.write("\n")

    def output_subsystem(self, f):
        """
        Output subsystem information based on roles.

        :param f: File object to print to
        :type f: file
        """
        # Get subsystem information
        roles_to_ss = PyFBA.parse.roles_to_subsystem(set(self.roles.keys()))
        f.write("function\tsubsystem\tsubcategory\tcategory\n")
        for role, info in roles_to_ss.items():
            for i in info:
                cat, subcat, ss = i
                f.write("{}\t{}\t{}\t{}\n".format(role, ss, subcat, cat))

    def reactions_to_subsystems(self, ss_file):
        """
        Output subsystem information for each reaction.
        First, roles that were used to build model are matched up with
        roles-subsystem relationship present in the RAST file. Then, if there
        are any other leftover roles, their subsystems are gathered from the
        PyFBA roles-subsystem file. Finally, using the reaction-role
        relationship used at time of model creation, reactions are paired up
        with subsystems.
        Returned dictionary has the reaction ID as the key and a list of
        4-tuples: (role, category, subcategory, subsystem).

        :param ss_file: Subsystem file from RAST
        :type ss_file: str
        :return: Mapping of reaction IDs to subsystems
        :rtype: dict
        """
        # Check if file exists
        if not os.path.isfile(ss_file):
            raise ValueError("Subsystem file not found")

        # Read in RAST subsystems file
        ss_file_data = dict()
        with open(ss_file, "r") as f:
            # Discard header
            f.readline()
            for l in f:
                l = l.rstrip("\n")
                cat, subcat, ss, role, features = l.split("\t")

                # Add subsystem info as a 3-tuple
                if role not in ss_file_data:
                    ss_file_data[role] = []
                to_add = (cat, subcat, ss)
                ss_file_data[role].append(to_add)

        # Match model roles to those from the file
        roles_to_ss = dict()
        model_roles = set(self.roles.keys())  # Temporary set
        for role in self.roles.keys():
            # Check if role was in file
            if role not in ss_file_data:
                continue

            if role not in roles_to_ss:
                roles_to_ss[role] = []
            # Add all of the file data
            roles_to_ss[role].extend(ss_file_data[role])
            # Remove from temporary set
            model_roles.remove(role)

        # Find subsystem information for roles not in file
        if len(model_roles) > 0:
            mr_to_ss = PyFBA.parse.roles_to_subsystem(model_roles)
            for role, info in mr_to_ss.items():
                roles_to_ss[role] = info

        # Get mapping from model reaction IDs to roles
        mReactions = {r: [] for r in self.reactions.keys()}
        for role, rxns in self.roles.items():
            for r in rxns:
                mReactions[r].append(role)

        # Initiate reactions
        rxn_to_ss = {rid: set() for rid in self.reactions}

        # Iterate through reactions and match to roles
        for rxn in rxn_to_ss:
            # Check if reaction has a role, it might not if it was gap-filled
            if rxn not in mReactions:
                data = ("Unknown", "Unknown", "Unknown", "Unknown")
                rxn_to_ss[rxn].add(data)
            else:
                for role in mReactions[rxn]:
                    # Check if role was found previously
                    if role not in roles_to_ss:
                        data = (role, "Unknown", "Unknown", "Unknown")
                        rxn_to_ss[rxn].add(data)
                    else:
                        for info in roles_to_ss[role]:
                            cat, subcat, ss = info
                            data = (role, cat, subcat, ss)
                            rxn_to_ss[rxn].add(data)

        return rxn_to_ss

    def run_fba(self, media_file, biomass_reaction=None, custom_bounds={},
                verbose=False):
        """
        Run FBA on model and return status, value, and growth.

        :param media_file: Media filepath
        :type media_file: str
        :param biomass_reaction: Given biomass Reaction object
        :type biomass_reaction: Reaction
        :param custom_bounds: A hash of Reactions (keys) and their bounds as a 2-tuple (values)
        :type custom_bounds: dict
        :return: FBA status, the biomass reaction flux, and
        :rtype: tuple of str, float, and bool
        """
        # Check if model has a biomass reaction if none was given
        if not biomass_reaction and not self.biomass_reaction:
            raise Exception("Model has no biomass reaction, "
                            "please supply one to run FBA")

        elif not biomass_reaction:
            biomass_reaction = self.biomass_reaction

        # Read in media file
        try:
            media = PyFBA.parse.read_media_file(media_file)
        except IOError as e:
            print(e)
            return None, None, None

        # Load ModelSEED database
        compounds, reactions, enzymes =\
            PyFBA.parse.model_seed.compounds_reactions_enzymes(
                self.organism_type)

        modelRxns = [rID for rID in self.reactions]
        modelRxns = set(modelRxns)

        status, value, growth = PyFBA.fba.run_fba(compounds,
                                                  reactions,
                                                  modelRxns,
                                                  media,
                                                  biomass_reaction,
                                                  custom_bounds=custom_bounds,
                                                  verbose=verbose)

        return status, value, growth

    def run_fba_pm_plate(self, pm_plate, mask=set(), biomass_reaction=None,
                         verbose=False):
        """
        Run FBA across a PM plate

        :param pm_plate: Name of the PM plate
        :type pm_plate: str
        :param mask: Set of wells to omit from the FBA runs
        :type mask: set
        :param biomass_reaction: Given biomass Reaction object
        :type biomass_reaction: Reaction
        :param verbose: Verbose output flag
        :type verbose: bool
        :return: Growth information for each media
        :rtype: dict of dict
        """
        results = {}
        # Set plate compositions
        pm_plates = PyFBA.parse.read_media.get_pm_plates()

        # Check if PM plate exist
        if pm_plate not in pm_plates:
            raise Exception("The PM plate " + pm_plate + " does not exist")

        # Iterate through media
        num_media = len(pm_plates[pm_plate])
        if verbose:
            print("Preparing to run FBA on a total of", num_media,
                  "media conditions")
        for midx, media in enumerate(pm_plates[pm_plate], start=1):
            # Check if media must be skipped
            if media in mask:
                if verbose:
                    print(media, "in mask set. Skipping.", file=sys.stderr)
                # Remove from mask set
                mask.remove(media)
                continue

            # Run the FBA
            if verbose:
                print("FBA running on", media,
                      "({} out of {})".format(midx, num_media),
                      file=sys.stderr)
            status, flux, growth = self.run_fba(media, biomass_reaction)
            if verbose:
                print(media, ": ", growth, " (", flux, ")",
                      sep="", file=sys.stderr)

            results[media] = {"status": status,
                              "biomass_flux": flux,
                              "growth": growth}

        if verbose:
            print("FBA on", pm_plate, "complete!", file=sys.stderr)

        return results

    def find_essential_reactions(self, media, verbose=False):
        """
        Find which reactions in the model are essential for growth on a given
        media. Essential reactions are those that, when removed from the model,
        result in no growth
        :param media: Media file
        :type media: str
        :param verbose: Verbose output flag
        :type verbose: bool
        :return: Essential reactions
        :rtype: set
        """
        essential = set()

        num = self.number_of_reactions()
        # Iterate through all reactions
        for i, rxn in enumerate(self.reactions, start=1):
            # Make copy of the model
            tmp_model = PyFBA.model.Model(self.id, self.name,
                                          self.organism_type)
            tmp_model.reactions = copy.deepcopy(self.reactions)
            tmp_model.compounds = copy.deepcopy(self.compounds)
            tmp_model.biomass_reaction = copy.copy(self.biomass_reaction)

            # Remove reaction from Model
            if verbose:
                eqn = self.reactions[rxn].equation
                print("(", i, "/", num, ") Removing ", rxn, ": ", eqn,
                      " Current count: ", len(essential),
                      file=sys.stderr, end="", sep="")
            del tmp_model.reactions[rxn]

            # Run FBA
            status, value, growth = tmp_model.run_fba(media)
            if not growth:
                essential.add(rxn)

        if verbose:
            print('Model contains', len(essential), 'essential reactions')

        return essential

    def find_essential_reactions_mt(self, media, n_threads=1, verbose=False):
        """
        Find which reactions in the model are essential for growth on a given
        media. Essential reactions are those that, when removed from the model,
        result in no growth. Function provides use of multithreading with the
        threading and queue modules.

        :param media: Media file
        :type media: str
        :param n_threads: Number of threads to use
        :type n_threads: int
        :param verbose: Verbose output flag
        :type verbose: bool
        :return: Essential reactions
        :rtype: set
        """
        ##### IMPORTANT #####
        # Currently, one thread is only allowed due to the setup of PyFBA.
        # A shared GLPK solver is used and may be compromised between threads
        # when multi-threading is implemented. This may be resolved in future
        # versions
        #####################
        import threading
        import queue
        thread_statement = "WARNING: Currently, only one thread is allowed " \
                           "due to the setup of PyFBA. A shared GLPK solver " \
                           "is used and may be compromised between threads " \
                           "when multi-threading is implemented. This may " \
                           "be resolved in future versions."
        if n_threads != 1:
            print(thread_statement, file=sys.stderr)
            n_threads = 1

        # Global variables
        essential = set()  # Set object of reactions to return
        rxn_q = queue.Queue()  # Main queue of reactions
        q_lock = threading.Lock()  # Queue lock for threads
        threads = []

        # Create Thread class
        class ER_Thread(threading.Thread):
            def __init__(self, thread_id, q):
                threading.Thread.__init__(self)
                self.thread_id = thread_id
                self.q = q

            def run(self):
                process_data(self.thread_id, self.q)

        # Create thread worker function
        def process_data(thread_id, q):
            """
            Perform infinite loop of FBA until reaction queue is completed

            :param thread_id: Thread ID
            :type thread_id: int
            :param q: Queue of reaction IDs
            :type q: Queue
            :return: None
            """
            while True:
                # Grab queue lock
                q_lock.acquire()

                if rxn_q.empty():
                    q_lock.release()
                    if verbose:
                        print("Thread", thread_id, "is closing",
                              file=sys.stderr)
                    break

                # When queue is not empty, grab next reaction ID
                r_id = q.get()
                q_lock.release()

                # Make copy of the model
                tmp_model = PyFBA.model.Model(self.id, self.name,
                                              self.organism_type)
                tmp_model.reactions = copy.deepcopy(self.reactions)
                tmp_model.compounds = copy.deepcopy(self.compounds)
                tmp_model.biomass_reaction = copy.copy(self.biomass_reaction)

                # Remove reaction from Model
                if verbose:
                    print("Thread", thread_id, ": Removing", r_id,
                          "from the model", file=sys.stderr)
                del tmp_model.reactions[r_id]

                # Run FBA
                status, value, growth = tmp_model.run_fba(media)
                if verbose:
                    ans = "not " if growth else ""
                    print("Thread ", thread_id, " : ", r_id, " is ", ans,
                          "essential", sep="", file=sys.stderr)
                if not growth:
                    essential.add(r_id)

        # Fill up queue
        for rxn in self.reactions:
            rxn_q.put(rxn)

        if verbose:
            print("Model contains", rxn_q.qsize(), "reactions to test",
                  file=sys.stderr)

        # Start threads
        if verbose:
            print("Starting up", n_threads, "threads", file=sys.stderr)
        for i in range(n_threads):
            t = ER_Thread(i, rxn_q)
            t.start()
            threads.append(t)

        # Synchronize threads -- wait for them to complete
        for t in threads:
            t.join()

        if verbose:
            print("Model contains", len(essential), "essential reactions")

        return essential

    def gapfill(self, media_file, cg_file, use_flux=False, verbose=0):
        """
        Gap-fill model on given media.

        :param media_file: Media filepath
        :type media_file: str
        :param cg_file: Close genomes roles filepath
        :type cg_file: str
        :param use_flux: Flag to remove candidate reactions carrying no flux\
        :type use_flux: bool
        :param verbose: Verbose output level
        :type verbose: int
        :rtype: bool
        """
        # Check if model needs any gap-filling
        status, value, growth = self.run_fba(media_file)

        # Check that FBA ran successfully
        if not status:
            return False
        elif growth:
            print("Initial FBA run results in growth; no need to gap-fill on given media")
            print("Biomass flux value:", value)
            sys.stdout.flush()
            # Add media to gap-filled media
            self.gapfilled_media.add(basename(media_file))
            return True


        # Read in media file
        try:
            media = PyFBA.parse.read_media_file(media_file)
        except IOError as e:
            print(e)
            return False

        # Create new model to run gap-filling with
        newModel = PyFBA.model.Model(self.id, self.name, self.organism_type)
        newModel.add_reactions(set(self.reactions.values()))
        newModel.set_biomass_reaction(self.biomass_reaction)

        newModelRxns = [rID for rID in self.reactions]
        newModelRxns = set(newModelRxns)
        original_reactions = copy.copy(newModelRxns)
        added_reactions = []

        if verbose >= 1:
            print("Current model contains", len(newModelRxns), "reactions",
                  file=sys.stderr)
            sys.stderr.flush()

        # Load ModelSEED database
        compounds, reactions, enzymes =\
            PyFBA.parse.model_seed.compounds_reactions_enzymes(
                self.organism_type)

        ########################################
        ## Media import reactions
        ########################################
        if verbose >= 1:
            print("Finding media import reactions", file=sys.stderr)
            sys.stderr.flush()
        gf_reactions = PyFBA.gapfill.suggest_from_media(compounds,
                                                        reactions,
                                                        newModelRxns,
                                                        media)
        added_reactions.append(("media", gf_reactions))
        newModelRxns.update(gf_reactions)
        rxns_for_new_model = set()
        for r in gf_reactions:
            rxns_for_new_model.add(reactions[r])
        newModel.add_reactions(rxns_for_new_model)
        if verbose >= 1:
            print("Found", len(gf_reactions), "reactions", file=sys.stderr)
            print("New total:", len(newModelRxns), "reactions", file=sys.stderr)
            sys.stderr.flush()

        if len(gf_reactions) > 0:
            # Run FBA
            status, value, growth = newModel.run_fba(media_file)
        if not growth:
            ####################################
            ## Essential reactions
            ####################################
            if verbose >= 1:
                print("Finding essential reactions", file=sys.stderr)
                sys.stderr.flush()
            gf_reactions = PyFBA.gapfill.suggest_essential_reactions()
            added_reactions.append(("essential", gf_reactions))
            newModelRxns.update(gf_reactions)
            rxns_for_new_model = set()
            for r in gf_reactions:
                rxns_for_new_model.add(reactions[r])
            newModel.add_reactions(rxns_for_new_model)
            if verbose >= 1:
                print("Found", len(gf_reactions), "reactions",
                      file=sys.stderr)
                print("New total:", len(newModelRxns), "reactions",
                      file=sys.stderr)
                sys.stderr.flush()

            if len(gf_reactions) > 0:
                # Run FBA
                status, value, growth = newModel.run_fba(media_file)
        if not growth:
            ####################################
            ## Close organism reactions
            ####################################
            if verbose >= 1:
                print("Finding close organism reactions", file=sys.stderr)
                sys.stderr.flush()
            gf_reactions =\
                    PyFBA.gapfill.suggest_from_roles(cg_file, reactions)
            added_reactions.append(("close genomes", gf_reactions))
            newModelRxns.update(gf_reactions)
            rxns_for_new_model = set()
            for r in gf_reactions:
                rxns_for_new_model.add(reactions[r])
            newModel.add_reactions(rxns_for_new_model)
            if verbose >= 1:
                print("Found", len(gf_reactions), "reactions",
                      file=sys.stderr)
                print("New total:", len(newModelRxns), "reactions",
                      file=sys.stderr)
                sys.stderr.flush()

            if len(gf_reactions) > 0:
                # Run FBA
                status, value, growth = newModel.run_fba(media_file)
        if not growth:
            ####################################
            ## Subsystem reactions
            ####################################
            if verbose >= 1:
                print("Finding subsystem reactions", file=sys.stderr)
                sys.stderr.flush()
            gf_reactions =\
                    PyFBA.gapfill.suggest_reactions_from_subsystems(reactions,
                                                                    newModelRxns,
                                                                    threshold=0.5)
            added_reactions.append(("subsystems", gf_reactions))
            newModelRxns.update(gf_reactions)
            rxns_for_new_model = set()
            for r in gf_reactions:
                rxns_for_new_model.add(reactions[r])
            newModel.add_reactions(rxns_for_new_model)
            if verbose >= 1:
                print("Found", len(gf_reactions), "reactions",
                      file=sys.stderr)
                print("New total:", len(newModelRxns), "reactions",
                      file=sys.stderr)
                sys.stderr.flush()

            if len(gf_reactions) > 0:
                # Run FBA
                status, value, growth = newModel.run_fba(media_file)
        if not growth:
            ####################################
            ## EC reactions
            ####################################
            if verbose >= 1:
                print("Finding EC reactions", file=sys.stderr)
                sys.stderr.flush()
            gf_reactions =\
                    PyFBA.gapfill.suggest_reactions_using_ec(newModel.roles.keys(),
                                                             reactions,
                                                             newModelRxns)
            added_reactions.append(("ec", gf_reactions))
            newModelRxns.update(gf_reactions)
            rxns_for_new_model = set()
            for r in gf_reactions:
                rxns_for_new_model.add(reactions[r])
            newModel.add_reactions(rxns_for_new_model)
            if verbose >= 1:
                print("Found", len(gf_reactions), "reactions",
                      file=sys.stderr)
                print("New total:", len(newModelRxns), "reactions",
                      file=sys.stderr)
                sys.stderr.flush()

            if len(gf_reactions) > 0:
                # Run FBA
                status, value, growth = newModel.run_fba(media_file)
        if not growth:
            ####################################
            ## Compound-probabilty reactions
            ####################################
            if verbose >= 1:
                print("Finding compound-probability reactions", file=sys.stderr)
                sys.stderr.flush()
            gf_reactions=\
                    PyFBA.gapfill.compound_probability(reactions,
                                                       newModelRxns,
                                                       cutoff=0,
                                                       rxn_with_proteins=True)
            added_reactions.append(("probable", gf_reactions))
            newModelRxns.update(gf_reactions)
            rxns_for_new_model = set()
            for r in gf_reactions:
                rxns_for_new_model.add(reactions[r])
            newModel.add_reactions(rxns_for_new_model)
            if verbose >= 1:
                print("Found", len(gf_reactions), "reactions",
                      file=sys.stderr)
                print("New total:", len(newModelRxns), "reactions",
                      file=sys.stderr)
                sys.stderr.flush()

            if len(gf_reactions) > 0:
                # Run FBA
                status, value, growth = newModel.run_fba(media_file)
        if not growth:
            ####################################
            ## Orphan-compound reactions
            ####################################
            if verbose >= 1:
                print("Finding orphan-compound reactions", file=sys.stderr)
                sys.stderr.flush()
            gf_reactions =\
                    PyFBA.gapfill.suggest_by_compound(compounds,
                                                      reactions,
                                                      newModelRxns,
                                                      max_reactions=1)
            added_reactions.append(("orphans", gf_reactions))
            newModelRxns.update(gf_reactions)
            rxns_for_new_model = set()
            for r in gf_reactions:
                rxns_for_new_model.add(reactions[r])
            newModel.add_reactions(rxns_for_new_model)
            if verbose >= 1:
                print("Found", len(gf_reactions), "reactions",
                      file=sys.stderr)
                print("New total:", len(newModelRxns), "reactions",
                      file=sys.stderr)
                sys.stderr.flush()

            if len(gf_reactions) > 0:
                # Run FBA
                status, value, growth = newModel.run_fba(media_file)
        ########################################
        ## Check if gap-filling was successful
        ########################################
        if not growth:
            print("Unable to gap-fill model", file=sys.stderr)
            sys.stderr.flush()
            return False
        elif verbose >= 1:
            print("Gap-fill was successful, now trimming model",
                  file=sys.stderr)
            sys.stderr.flush()

        ########################################
        ## Trimming the model
        ########################################
        if use_flux:
            # Get fluxes from gap-filled reactions
            # Keep those without a flux of zero
            rxnfluxes = PyFBA.model.model_reaction_fluxes(newModel,
                                                          media_file)
            numRemoved = 0
            tmp_added_reactions = []
            for how, gfrxns in added_reactions:
                tmp_gfrxns = set()
                for gfr in gfrxns:
                    if float(rxnfluxes[gfr]) == 0.0:
                        numRemoved += 1
                    else:
                        tmp_gfrxns.add(gfr)
                tmp_added_reactions.append((how, tmp_gfrxns))
            added_reactions = tmp_added_reactions

            if verbose >= 1:
                print("Removed", numRemoved, "reactions based on flux value",
                    file=sys.stderr)

        required_rxns = set()
        gapfilled_keep = set()
        # Begin loop through all gap-filled reactions
        while added_reactions:
            ori = copy.copy(original_reactions)
            ori.update(required_rxns)

            # Test next set of gap-filled reactions
            # Each set is based on a method described above
            how, new = added_reactions.pop()

            if verbose >= 1:
                print("Trimming", how, "group of reactions",
                      file=sys.stderr)
                sys.stderr.flush()

            # Get all the other gap-filled reactions we need to add
            for gf_tple in added_reactions:
                ori.update(gf_tple[1])

            # Use minimization function to determine the minimal
            # set of gap-filled reactions from the current method
            verb = verbose == 2
            minimized_set =\
                    PyFBA.gapfill.minimize_additional_reactions(ori,
                                                                new,
                                                                compounds,
                                                                reactions,
                                                                media,
                                                                newModel.biomass_reaction,
                                                                verbose=verb)
            # Record the method used to determine
            # how the reaction was gap-filled
            for new_rxn in minimized_set:
                reactions[new_rxn].is_gapfilled = True
                reactions[new_rxn].gapfill_method = how
            required_rxns.update(minimized_set)
            gapfilled_keep.update(minimized_set)
        # End trimming

        if verbose >= 1:
            print("Trimming complete.\nTotal gap-filled reactions:",
                  len(required_rxns), file=sys.stderr)
            sys.stderr.flush()

        # Record reactions and roles for each gap-filled reaction
        add_to_model_rxns = set()
        add_to_model_roles = {}
        gf_reactions = PyFBA.filters.reactions_to_roles(gapfilled_keep, verb)
        for rxn in gapfilled_keep:
            if rxn in original_reactions:
                continue
            # Add to model gf_reactions dictionary
            self.gf_reactions[rxn] = (reactions[rxn].gapfill_method,
                                      basename(media_file))
            add_to_model_rxns.add(reactions[rxn])
            try:
                for rl in gf_reactions[rxn]:
                    if rl not in add_to_model_roles:
                        add_to_model_roles[rl] = set()
                    add_to_model_roles[rl].add(rxn)
            except KeyError:
                pass

        # Add to this model
        self.add_reactions(add_to_model_rxns)
        self.add_roles(add_to_model_roles)

        # Add media to gap-filled media
        self.gapfilled_media.add(basename(media_file))

        # Run FBA
        status, value, growth = self.run_fba(media_file)
        if not growth:
            print("Failed final FBA check!", file=sys.stderr)
            return False
        print("The biomass reaction has a flux of", value)

        return growth
