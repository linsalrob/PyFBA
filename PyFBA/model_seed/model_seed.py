"""

"""

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

    def reset(self):
        self.compounds = None
        self.reactions = {}
        self.enzymes = None
