from __future__ import print_function
import PyFBA

class Model:
    """
    A model is a structured collection of reactions, compounds, and functional roles.

    A model will provide methods for extracting its information and creating SBML-formatted files.

    :ivar id: Model ID
    :ivar name: Model name
    :ivar roles: A dictionary of roles as the key and sets of reaction IDs as the value
    :ivar reactions: A set of reaction objects contained in this model
    :ivar compounds: A set of compound objects contained in this model
    :ivar gapfilled_media: A set of media names this model has been gap-filled with
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
        self.reactions = set()
        self.compounds = set()
        self.gapfilled_media = set()
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
        Add a role to the model.

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
        Add a reaction to the model.

        :param rxns: A Reaction object
        :type rxns: set
        """
        if isinstance(rxns, set):
            self.reactions.update(rxns)
            # Add compounds from new reactions
            for r in rxns:
                self.compounds.update(r.all_compounds())
        else:
            raise TypeError("You need to add a set of reactions to a model")


    def has_reaction(self, rxn):
        """
        Check if model contains reaction.

        :param rxn: A Reaction object
        :type rxn: Reaction
        :rtype: bool
        """
        return rxn in self.reactions


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
            return

        elif not biomass_reaction:
            biomass_reaction = self.biomass_reaction

        # Read in media file
        try:
            media = PyFBA.parse.read_media_file(media_file)
        except IOError as e:
            print(e)
            return

        # Load ModelSEED database
        compounds, reactions, enzymes =\
            PyFBA.parse.model_seed.compounds_reactions_enzymes(
                self.organism_type)

        modelRxns = [r.name for r in self.reactions]
        modelRxns = set(modelRxns)

        status, value, growth = PyFBA.fba.run_fba(compounds,
                                                  reactions,
                                                  modelRxns,
                                                  media,
                                                  biomass_reaction)

        return (status, value, growth)
