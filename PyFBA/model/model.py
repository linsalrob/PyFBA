from __future__ import print_function

from typing import Set, Any

import PyFBA
import copy
import sys
from os.path import basename


class Model:
    """
    A model is a structured collection of reactions, compounds, and functional roles.

    A model will provide methods for extracting its information and creating SBML-formatted files.

    :ivar id: Model ID
    :ivar name: Model name
    :ivar roles: A dictionary of roles as the key and sets of reaction IDs as the value
    :ivar reactions: A dictionary of reaction IDs as the key and Reaction objects as the value
    :ivar compounds: A dictionary of compound IDs as the key and Compound objects as the value
    :ivar gapfilled_media: A set of media names this model has been gap-filled with
    :ivar gf_reactions: A set of gap-filled reaction IDs
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
        self.roles = {}
        self.reactions = {}
        self.compounds = {}
        self.gapfilled_media = set()
        self.gf_reactions = set()
        self.biomass_reaction = None
        self.organism_type = organism_type

    def __str__(self):
        """
        The to string function.

        :rtype: str
        """
        return self.name + " (id: " + self.id + ")"
    
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
                if r.id not in self.reactions:
                    self.reactions[r.id] = r
                    # Add compounds from new reactions
                    for c in r.all_compounds():
                        if c.id not in self.compounds:
                            self.compounds[c.id] = c
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
        return cpd.id in self.compounds

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
        m_reactions = {r: [] for r in self.reactions.keys()}
        for role, rxns in self.roles.items():
            for r in rxns:
                m_reactions[r].append(role)
        # Print header
        f.write("reaction\tfunction\tequation\tgapfilled\n")

        # Print reactions from model
        for r, roles in m_reactions.items():
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

    def run_fba(self, media_file, biomass_reaction=None):
        """
        Run FBA on model and return status, value, and growth.

        :param media_file: Media filepath
        :type media_file: str
        :param biomass_reaction: Given biomass Reaction object
        :type biomass_reaction: Reaction
        :rtype: tuple
        """
        # Check if model has a biomass reaction if none was given
        if not biomass_reaction and not self.biomass_reaction:
            raise Exception("Model has no biomass reaction, please supply one to run FBA")
        elif not biomass_reaction:
            biomass_reaction = self.biomass_reaction

        # Read in media file
        try:
            media = PyFBA.parse.read_media_file(media_file)
        except IOError as e:
            print(e)
            return None, None, None

        # Load ModelSEED database
        modeldata = PyFBA.parse.model_seed.parse_model_seed_data(self.organism_type, verbose=False)

        model_rxns = [rID for rID in self.reactions]
        model_rxns = set(model_rxns)

        status, value, growth = PyFBA.fba.run_fba(modeldata,
                                                  model_rxns,
                                                  media,
                                                  biomass_reaction)

        return status, value, growth

    def gapfill(self, media_file, cg_file, use_flux=False, verbose=0):
        """
        Gap-fill model on given media.

        :param use_flux: Use the flux
        :param media_file: Media filepath
        :type media_file: str
        :param cg_file: Close genomes roles filepath
        :type cg_file: str
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
        new_model = PyFBA.model.Model(self.id, self.name, self.organism_type)
        new_model.add_reactions(set(self.reactions.values()))
        new_model.set_biomass_reaction(self.biomass_reaction)

        new_model_rxns = [rID for rID in self.reactions]
        new_model_rxns = set(new_model_rxns)
        original_reactions = copy.copy(new_model_rxns)
        added_reactions = []

        if verbose >= 1:
            print("Current model contains", len(new_model_rxns), "reactions",
                  file=sys.stderr)
            sys.stderr.flush()

        # Load ModelSEED database
        modeldata = PyFBA.parse.model_seed.parse_model_seed_data(self.organism_type, verbose=False)

        ########################################
        # Media import reactions
        ########################################
        if verbose >= 1:
            print("Finding media import reactions", file=sys.stderr)
            sys.stderr.flush()
        gf_reactions = PyFBA.gapfill.suggest_from_media(modeldata,
                                                        new_model_rxns,
                                                        media)
        added_reactions.append(("media", gf_reactions))
        new_model_rxns.update(gf_reactions)
        rxns_for_new_model = set()
        for r in gf_reactions:
            rxns_for_new_model.add(modeldata.reactions[r])
        new_model.add_reactions(rxns_for_new_model)
        if verbose >= 1:
            print("Found", len(gf_reactions), "reactions", file=sys.stderr)
            print("New total:", len(new_model_rxns), "reactions", file=sys.stderr)
            sys.stderr.flush()

        if len(gf_reactions) > 0:
            # Run FBA
            status, value, growth = new_model.run_fba(media_file)
        if not growth:
            ####################################
            # Essential reactions
            ####################################
            if verbose >= 1:
                print("Finding essential reactions", file=sys.stderr)
                sys.stderr.flush()
            gf_reactions = PyFBA.gapfill.suggest_essential_reactions()
            added_reactions.append(("essential", gf_reactions))
            new_model_rxns.update(gf_reactions)
            rxns_for_new_model = set()
            for r in gf_reactions:
                rxns_for_new_model.add(modeldata.reactions[r])
            new_model.add_reactions(rxns_for_new_model)
            if verbose >= 1:
                print("Found", len(gf_reactions), "reactions",
                      file=sys.stderr)
                print("New total:", len(new_model_rxns), "reactions",
                      file=sys.stderr)
                sys.stderr.flush()

            if len(gf_reactions) > 0:
                # Run FBA
                status, value, growth = new_model.run_fba(media_file)
        if not growth:
            ####################################
            # Close organism reactions
            ####################################
            if verbose >= 1:
                print("Finding close organism reactions", file=sys.stderr)
                sys.stderr.flush()
            gf_reactions = PyFBA.gapfill.suggest_from_roles(cg_file, modeldata.reactions)
            added_reactions.append(("close genomes", gf_reactions))
            new_model_rxns.update(gf_reactions)
            rxns_for_new_model = set()
            for r in gf_reactions:
                rxns_for_new_model.add(modeldata.reactions[r])
            new_model.add_reactions(rxns_for_new_model)
            if verbose >= 1:
                print("Found", len(gf_reactions), "reactions",
                      file=sys.stderr)
                print("New total:", len(new_model_rxns), "reactions",
                      file=sys.stderr)
                sys.stderr.flush()

            if len(gf_reactions) > 0:
                # Run FBA
                status, value, growth = new_model.run_fba(media_file)
        if not growth:
            ####################################
            # Subsystem reactions
            ####################################
            if verbose >= 1:
                print("Finding subsystem reactions", file=sys.stderr)
                sys.stderr.flush()
            gf_reactions = PyFBA.gapfill.suggest_reactions_from_subsystems(modeldata.reactions, new_model_rxns, 
                                                                           threshold=0.5)
            added_reactions.append(("subsystems", gf_reactions))
            new_model_rxns.update(gf_reactions)
            rxns_for_new_model = set()
            for r in gf_reactions:
                rxns_for_new_model.add(modeldata.reactions[r])
            new_model.add_reactions(rxns_for_new_model)
            if verbose >= 1:
                print("Found", len(gf_reactions), "reactions",
                      file=sys.stderr)
                print("New total:", len(new_model_rxns), "reactions",
                      file=sys.stderr)
                sys.stderr.flush()

            if len(gf_reactions) > 0:
                # Run FBA
                status, value, growth = new_model.run_fba(media_file)
        if not growth:
            ####################################
            # EC reactions
            ####################################
            if verbose >= 1:
                print("Finding EC reactions", file=sys.stderr)
                sys.stderr.flush()
            gf_reactions = PyFBA.gapfill.suggest_reactions_using_ec(set(new_model.roles.keys()), modeldata.reactions,
                                                                    new_model_rxns)
            added_reactions.append(("ec", gf_reactions))
            new_model_rxns.update(gf_reactions)
            rxns_for_new_model = set()
            for r in gf_reactions:
                rxns_for_new_model.add(modeldata.reactions[r])
            new_model.add_reactions(rxns_for_new_model)
            if verbose >= 1:
                print("Found", len(gf_reactions), "reactions",
                      file=sys.stderr)
                print("New total:", len(new_model_rxns), "reactions",
                      file=sys.stderr)
                sys.stderr.flush()

            if len(gf_reactions) > 0:
                # Run FBA
                status, value, growth = new_model.run_fba(media_file)
        if not growth:
            ####################################
            # Compound-probabilty reactions
            ####################################
            if verbose >= 1:
                print("Finding compound-probability reactions", file=sys.stderr)
                sys.stderr.flush()
            gf_reactions = PyFBA.gapfill.compound_probability(modeldata.reactions, new_model_rxns, cutoff=0,
                                                              rxn_with_proteins=True)
            added_reactions.append(("probable", gf_reactions))
            new_model_rxns.update(gf_reactions)
            rxns_for_new_model = set()
            for r in gf_reactions:
                rxns_for_new_model.add(modeldata.reactions[r])
            new_model.add_reactions(rxns_for_new_model)
            if verbose >= 1:
                print("Found", len(gf_reactions), "reactions",
                      file=sys.stderr)
                print("New total:", len(new_model_rxns), "reactions",
                      file=sys.stderr)
                sys.stderr.flush()

            if len(gf_reactions) > 0:
                # Run FBA
                status, value, growth = new_model.run_fba(media_file)
                
        if not growth:
            ####################################
            # Orphan-compound reactions
            ####################################
            if verbose >= 1:
                print("Finding orphan-compound reactions", file=sys.stderr)
                sys.stderr.flush()
            gf_reactions = PyFBA.gapfill.suggest_by_compound(modeldata, new_model_rxns, max_reactions=1)
            added_reactions.append(("orphans", gf_reactions))
            new_model_rxns.update(gf_reactions)
            rxns_for_new_model = set()
            for r in gf_reactions:
                rxns_for_new_model.add(modeldata.reactions[r])
            new_model.add_reactions(rxns_for_new_model)
            if verbose >= 1:
                print("Found", len(gf_reactions), "reactions",
                      file=sys.stderr)
                print("New total:", len(new_model_rxns), "reactions",
                      file=sys.stderr)
                sys.stderr.flush()

            if len(gf_reactions) > 0:
                # Run FBA
                status, value, growth = new_model.run_fba(media_file)
       
        ########################################
        # Check if gap-filling was successful
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
        # Trimming the model
        ########################################
        if use_flux:
            # Get fluxes from gap-filled reactions
            # Keep those without a flux of zero
            rxnfluxes = PyFBA.model.model_reaction_fluxes(new_model,
                                                          media_file)
            num_removed = 0
            tmp_added_reactions = []
            for how, gfrxns in added_reactions:
                tmp_gfrxns = set()
                for gfr in gfrxns:
                    if float(rxnfluxes[gfr]) == 0.0:
                        num_removed += 1
                    else:
                        tmp_gfrxns.add(gfr)
                tmp_added_reactions.append((how, tmp_gfrxns))
            added_reactions = tmp_added_reactions

            if verbose >= 1:
                print("Removed", num_removed, "reactions based on flux value", file=sys.stderr)

        required_rxns = set()
        gapfilled_keep = set()
        # Begin loop through all gap-filled reactions
        verb = verbose == 2
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
            minimized_set = PyFBA.gapfill.minimize_additional_reactions(ori, new, modeldata, media,
                                                                        new_model.biomass_reaction, verbose=verb)
            # Record the method used to determine
            # how the reaction was gap-filled
            for new_rxn in minimized_set:
                modeldata.reactions[new_rxn].is_gapfilled = True
                modeldata.reactions[new_rxn].gapfill_method = how
            required_rxns.update(minimized_set)
            gapfilled_keep.update(minimized_set)
        # End trimming

        if verbose >= 1:
            print("Trimming complete.\nTotal gap-filled reactions:",
                  len(required_rxns), file=sys.stderr)
            sys.stderr.flush()

        # Record reactions and roles for each gap-filled reaction
        add_to_model_rxns: Set[PyFBA.metabolism.Reaction] = set()
        add_to_model_roles = {}
        gf_reactions = PyFBA.filters.reactions_to_roles(gapfilled_keep, verbose=verb)
        for rxn in gapfilled_keep:
            if rxn in original_reactions:
                continue
            self.gf_reactions.add(rxn)  # Add to model gf_reactions set
            add_to_model_rxns.add(modeldata.reactions[rxn])
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
