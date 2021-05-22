"""

"""
import PyFBA

class ModelSeed:
    """
     A class to hold model seed objects so that we only need to parse them once.

     This holds enzymes, reactions, and compounds. The variables are the same data structures
     returned by the parser (see parse/model_seed.py).

     Other variables associated with the Compound class:
     :ivar compounds: a dict of compound id -> compound objects
     :ivar reactions: a dict of organism type -> dict(reaction id -> reaction objects).
     :ivar enzymes:a dict of enzyme id -> enzyme objects

     """

    def __init__(self, compounds=None, reactions=None, enzymes=None):
        self.compounds = compounds
        if reactions:
            self.reactions = reactions
        else:
            self.reactions = {}
        self.enzymes = enzymes
        self.compounds_by_name = {}
        self.last_compound_by_name_sz = 0

    def reset(self):
        self.compounds = None
        self.reactions = {}
        self.enzymes = None

    def get_compound_by_name(self, name) -> PyFBA.metabolism.Compound:
        """
        Retrieve a compound by its name. We use self.last_compound_by_name_sz to see if compounds has changed
        :param name: The name to look through
        :return: the compound if found or None
        """

        if self.last_compound_by_name_sz != len(self.compounds):
            for c in self.compounds:
                self.compounds_by_name[self.compounds[c].name] = c
            self.last_compound_by_name_sz = len(self.compounds)

        if name in self.compounds_by_name:
            return self.compounds_by_name[name]
        return None
