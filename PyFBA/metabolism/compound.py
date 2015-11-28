
COMMON_REACTION_LIMIT = 5


class Compound:
    """
    A compound is the essential metabolic compound that is involved in a reaction.

    A compound has at the very minimum a name and a location. The location is typically one of:
       * e: extracellular
       * c: cytoplasmic
       * h: chloroplast
       * p: periplasm

    Other variables associated with the Compound class:
    :ivar name: the name of the  compound
    :ivar location: the location of the compound.
    :ivar reactions: a set of reaction objects that this compound is connected to
    :ivar model_seed_id: the compound id from the model seed.
    :ivar abbreviation: a short name for the compound
    :ivar formula: the compounds formula
    :ivar mw: the molecular weight of the compound
    :ivar common: Boolean: this is a common compound. This means the coompound is in > COMMON_REACTION_LIMIT reactions
    :ivar charge: the charge associated with the compound
    :ivar uptake_secretion: The compound is involved in uptake from the media or secretion back to the media

    """

    def __init__(self, name, location):
        """
        Initiate the object

        :param name: The name of the compound
        :type name: str
        :param location: The location of the compound
        :type location: str
        :return:
        :rtype:
        """
        self.name = name
        self.location = location
        self.reactions = set()
        self.model_seed_id = name
        self.alternate_seed_ids = set()
        self.abbreviation = None
        self.formula = None
        self.mw = 0
        self.common = False
        self.charge = 0
        self.uptake_secretion = False

    def __eq__(self, other):
        """
        Two compounds are equal if they have the same name and the same location

        :param other: The other compound
        :type other: Compound
        :return: If they are equal
        :rtype: bool
        """
        if isinstance(other, Compound):
            return (self.name, self.location) == (other.name, other.location)
        else:
            return NotImplemented
    
    def __cmp__(self, other):
        """
        Compare whether two things are the same.

        :param other: The other compound
        :type other: Compound
        :return: An int, zero if they are the same
        :rtype: int
        """
        if isinstance(other, Compound):
            if __eq__(other):
                return 0
            else:
                return 1
        else:
            return NotImplemented

    def __ne__(self, other):
        """
        Are these not equal?

        :param other: The other compound
        :type other: Compound
        :return: If they are not equal
        :rtype: bool
        """
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __hash__(self):
        """
        The hash function is based on the name of the compound.

        :rtype: int
        """
        return hash((self.name, self.location))

    def __str__(self):
        """
        The to string function.
        :rtype: str
        """
        return self.name + " (location: " + self.location + ")"


    def add_reactions(self, rxns):
        """
        Add a reaction that this compound is involved in. You  can add a set of reactions. See the note above about the
        number of reactions.

        :param rxns: A set of reactions
        :type rxns: set
        """
        if isinstance(rxns, set):
            self.reactions.update(rxns)
        else:
            raise TypeError("You need to add a set of reactions to a compound")


    def has_reaction(self, rxn):
        """
        Is this compound involved in this reaction?

        :param rxn: A Reaction object
        :type rxn: Reaction
        :return: Whether the reaction is present
        :rtype: bool
        """
        return rxn in self.reactions


    def number_of_reactions(self):
        """
        How many reactions is this compound involved in?

        :rtype: int
        """
        return len(self.reactions)


    def all_reactions(self):
        """
        Return a set of all the reactions that this compound is involved in

        :rtype: int
        """
        return self.reactions


    def is_common(self, rct_limit=COMMON_REACTION_LIMIT):
        """
        Is this a common compound? This requires that you have
        added reactions to this compound.
        
        You can either specify the number of reactions or use our
        default that is currently 50.

        :param rct_limit: The limit for a compound to be considered common
        :type rct_limit: int
        :return: Whether this is a common reaction
        :rtype: bool
        """
        if self.number_of_reactions() > rct_limit:
            self.common = True
        else:
            self.common = False
        return self.common

    def calculate_molecular_weight(self):
        """
        Calculate and return the molecular weight of this compound

        :return: The molecular weight
        :rtype: float
        """

        raise NotImplemented("Sorry. Calculate molecular weight has not yet been implemented.")