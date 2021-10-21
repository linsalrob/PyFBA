import sys

import PyFBA.metabolism


class Reaction:
    """
    A reaction is the central concept of metabolism and is the conversion of substrates to products.

    The reaction describes what we know. At a bare minimum we need a a name for the reaction. The name can either be the
    reaction id (e.g. modelSEED or KEGG id), or another name for this reaction.

    A reaction is an object that describes how to get from one compound to another. We need to know what the compound(s)
    on the left of the equation are, what the compounds on the right of the reaction are, and the probability that the
    reaction proceeds in either direction. If the reaction is truly reversible the probability can be 1 in both cases.
    If it is unidirectional the probability can be 0 in one direction.

    The likelihood that a reaction completes will be some product of its delta G and its p. We could also do something
    simpler, e.g. if there is a -ve delta G (favorable reaction) we can increase p and if there is a +ve delta G
    (unfavorable reaction) we can decrease p.

    The direction and reversible is the direction that the equation can run.

    Acceptable values are:
    ======   ===========================
    Value    Meaning
    ======   ===========================
    None     We don't know the direction
    >        Left to right
    <        Right to left
    =        Bidirectional
    ======   ===========================

    :ivar rctn_id: The reaction ID
    :ivar readable_name: The name of the reaction
    :ivar description: A description of the reaction
    :ivar equation: The reaction equation
    :ivar direction: The direction of the reaction (<, =, >, or ?)
    :ivar gfdirection: The possible gapfilled direction
    :ivar ntdirection: The non-template direction (before correcting for templates)
    :ivar left_compounds: A set of CompoundWithLocations on the left side of the reaction
    :ivar left_abundance: A dict of the CompoundWithLocations on the left and their abundance
    :ivar right_compounds: The set of CompoundWithLocations on the right side of the equation
    :ivar right_abundance: A dict of the CompoundWithLocations on the right and their abundance
    :ivar lower_bound: The lower bound for the reaction
    :ivar upper_bound: The upper bound for the reaction
    :ivar pLR: The probability the reaction proceeds left to right
    :ivar pRL: The probability the reaction proceeds right to left
    :ivar enzymes: The enzyme complex IDs involved in the reaction
    :ivar pegs: The protein-encoding genes involved in the reaction
    :ivar deltaG: The delta G
    :ivar deltaG_error: The error in the delta G
    :ivar inp: Whether the reaction is an input reaction
    :ivar outp: Whether the reaction is an output reaction
    :ivar is_transport: Whether the reaction is a transport reaction (imports or exports something)
    :ivar ran: Boolean to note whether the reaction ran
    :ivar is_biomass_reaction: Boolean to note whether this is a biomass reaction
    :ivar biomass_direction: If it is a biomass reaction, what is the direction
    :ivar is_gapfilled: Boolean to note whether the reaction was gapfilled
    :ivar gapfill_method: If the reaction was gapfilled, how was it gapfilled
    :ivar is_uptake_secretion: Is the reaction involved in uptake of compounds or secretion of compounds.

    """

    def __init__(self, rctn_id, readable_name=None, description=None, equation=None, direction=None):
        """
        Instantiate a reaction

        :param rctn_id: the reaction id
        :param readable_name: a human readable name. This was refactored from name to make it more unique
        :param description: a description of the reaction
        :param equation: the equation for the reaction
        :param direction: the direction of the reaction
        """
        self.id = rctn_id
        self.model_seed_id = rctn_id
        self.readable_name = readable_name
        self.description = description
        self.equation = equation
        self.direction = direction
        self.gfdirection = direction # the gap filling direction
        self.ntdirection = direction # the non-template driven direction
        self.left_compounds = set()  # type: set[PyFBA.metabolism.CompoundWithLocation]
        self.left_abundance = {}
        self.right_compounds = set()  # type: set[PyFBA.metabolism.CompoundWithLocation]
        self.right_abundance = {}
        self.lower_bound = None
        self.upper_bound = None
        self.pLR = 0
        self.pRL = 0
        self.enzymes = set()
        self.ec_numbers = []
        self.pegs = set()
        self.deltaG_error = 0
        self.deltaG = 0
        self.inp = False
        self.outp = False
        self.is_transport = False
        self.ran = False
        self.is_biomass_reaction = False
        self.biomass_direction = False
        self.is_gapfilled = False
        self.gapfill_method = ""
        self.is_uptake_secretion = False
        self.aliases = []

    def __eq__(self, other):
        """
        Two reactions are the same if they have the same left and
        right products, but not necessarily the same names or reactions.
        Note that we don't care whether the left and right (the 
        directionality) is the same in our two comparisons
        :param other: The other reaction
        :type other: Reaction
        :return: Boolean
        :rtype: bool
        """

        if isinstance(other, Reaction):
            return (self.id == other.id or
                    (self.left_compounds, self.right_compounds) ==
                    (other.left_compounds, other.right_compounds) or
                    (self.left_compounds, self.right_compounds) ==
                    (other.right_compounds, other.left_compounds)
                    )
        else:
            raise NotImplementedError(f"Comparing Reaction with {type(other)} is not implemented")

    def __cmp__(self, other):
        """
        Compare whether two things are the same
        :param other: The other reaction
        :type other: Reaction
        :return: an integer, zero if they are the same
        :rtype: int
        """

        if isinstance(other, Reaction):
            if __eq__(other):
                return 0
            else:
                return 1
        else:
            raise NotImplementedError(f"Comparing Reaction with {type(other)} is not implemented")

    def __ne__(self, other):
        """
        Are these not equal?
        :param other: The other reaction
        :type other: Reaction
        :return: Boolean
        :rtype: bool
        """

        try:
            result = self.__eq__(other)
        except NotImplementedError:
            return True
        return not result

    def __hash__(self):
        """
        The hash function is based on the name of the reaction.
        :rtype: int
        """
        return hash((self.id, self.readable_name))

    def __str__(self):
        """
        The string version of the reaction.
        :rtype: str
        """
        if self.readable_name:
            return f"{self.id}: {self.readable_name}"
        else:
            return f"{self.id}: {self.equation}"

    """
        Since we have complex data structures, we can't just pickle them and unpickle them with aplomb!
        In fact, this is affecting deep/shallow copy, and we need to ensure that we use copy.deepcopy() 
        at all times, otherwise the data structures are not copied correctly.
        
        These two methods correctly allow us to pickle the data structures. Note that we have
        CompoundWithLocation objects, and we need both the object and its abundance to correctly create the pickle.
    """

    def __getstate__(self):
        """
        The state that the object is saved or copied as. We override the left/right compounds and abundances
        with simple arrays of data. This is lossy - we are losing the connections between compounds and data
        and we probably need to reconstruct that after pickling/unpickling the reactions.
        :return:
        """
        state = self.__dict__.copy()
        state['left_compounds'] = []
        state['right_compounds'] = []
        state['left_abundance'] = {}
        state['right_abundance'] = {}
        for l in self.left_compounds:
            state['left_compounds'].append([l.id, l.name, l.location])
            state['left_abundance'][f"{l.id} :: {l.name} :: {l.location}"] = self.left_abundance[l]
        for r in self.right_compounds:
            state['right_compounds'].append([r.id, r.name, r.location])
            state['right_abundance'][f"{r.id} :: {r.name} :: {r.location}"] = self.right_abundance[r]
        return state


    def __setstate__(self, state):
        """
        Create a new reaction from a saved state. This is from __getstate__ eg. when pickled.
        :param state: the state that was saved.
        :return:
        """
        left = set()
        right = set()
        left_abundance = {}
        right_abundance = {}
        for l in state['left_compounds']:
            c = PyFBA.metabolism.CompoundWithLocation(id=l[0], name=l[1], location=l[2])
            left.add(c)
            left_abundance[c] = state['left_abundance'][f"{l[0]} :: {l[1]} :: {l[2]}"]
        state['left_compounds'] = left
        state['left_abundance'] = left_abundance
        for r in state['right_compounds']:
            c = PyFBA.metabolism.CompoundWithLocation(id=r[0], name=r[1], location=r[2])
            right.add(c)
            right_abundance[c] = state['right_abundance'][f"{r[0]} :: {r[1]} :: {r[2]}"]
        state['right_compounds'] = right
        state['right_abundance'] = right_abundance
        self.__dict__.update(state)



    def set_direction(self, direction):
        """
        Set the direction of the reaction.

        :param direction: The direction of the reaction
        :type direction: str
        :rtype: str
        :return: The current direction
        """

        allowable_directions = {'>', '<', '=', None}
        if direction in allowable_directions:
            self.direction = direction
            if not self.gfdirection:
                self.gfdirection = direction
        else:
            sys.stderr.write("Direction: " + str(direction) + " is not a permitted direction. Ignored\n")
            self.direction = None

        return self.direction

    def add_left_compounds(self, cmpds):
        """
        The compounds on the left are a set of compounds that the reaction typically uses as substrates.
        :param cmpds: The compounds that should be added
        :type cmpds: set[PyFBA.metabolism.CompoundWithLocation]
        """

        if isinstance(cmpds, set):
            # choose one element. next(iter(cmpds)) does not remove the element
            if not isinstance(next(iter(cmpds)), PyFBA.metabolism.CompoundWithLocation):
                raise TypeError(f"Starting with v.2 reactions need PyFBA.metabolism.CompoundWithLocation objects not {type(next(iter(cmpds)))}")
            self.left_compounds.update(cmpds)
        elif isinstance(cmpds, PyFBA.metabolism.CompoundWithLocation):
            # add a single compound
            self.left_compounds.add(cmpds)
        else:
            raise TypeError("Compounds must be a set of CompoundWithLocation")

    def set_left_compound_abundance(self, cmpd, abundance):
        """
        Set the abundance of a compound on the left side of the equation.

        :param cmpd: The compound to set the abundance for
        :type cmpd: PyFBA.metabolism.CompoundWithLocation
        :param abundance: The amount of that abundance
        :type abundance: float | int
        """

        if cmpd not in self.left_compounds:
            raise KeyError(f"{cmpd} is not in left compounds. Please add it before trying to set the abundance")
        if isinstance(abundance, float):
            self.left_abundance[cmpd] = abundance
        elif isinstance(abundance, int):
            self.left_abundance[cmpd] = float(abundance)
        else:
            raise TypeError("Abundance must be an int or a float")

    def get_left_compound_abundance(self, cmpd):
        """
        Get the abundance of the compound on the left side of the equation.

        :param cmpd: The compound to get the abundance of
        :type cmpd: PyFBA.metabolism.CompoundWithLocation
        :return: The compounds abundance
        :rtype: float
        """

        if cmpd in self.left_abundance:
            return self.left_abundance[cmpd]
        else:
            raise KeyError(f"In the reaction {self.readable_name} (reaction id: {self.id}), you do not have" +
                           f" {cmpd} on the left hand side of the equation: {self.equation}")

    def number_of_left_compounds(self):
        """
        The number of compounds on the left side of the equation.

        :rtype: int
        """
        return len(self.left_compounds)

    def add_right_compounds(self, cmpds):
        """
        The compounds on the right are a set of compounds that the reaction typically uses as substrates.

        :param cmpds: The compounds that should be added
        :type cmpds: set[PyFBA.metabolism.CompoundWithLocation]
        """
        if isinstance(cmpds, set):
            # choose one element. next(iter(cmpds)) does not remove the element
            if not isinstance(next(iter(cmpds)), PyFBA.metabolism.CompoundWithLocation):
                raise TypeError("Starting with v.2 reactions need PyFBA.metabolism.CompoundWithLocation objects")
            self.right_compounds.update(cmpds)
        elif isinstance(cmpds, PyFBA.metabolism.CompoundWithLocation):
            # add a single compound
            self.right_compounds.add(cmpds)
        else:
            raise TypeError("Compounds must be a set of CompoundWithLocation")

    def set_right_compound_abundance(self, cmpd, abundance):
        """
        Set the abundance of a compound on the right side of the equation

        :param cmpd: The compound to set the abundance for
        :type cmpd: PyFBA.metabolism.CompoundWithLocation
        :param abundance: The amount of that abundance
        :type abundance: float | int
        """
        if cmpd not in self.right_compounds:
            raise KeyError(f"{cmpd} is not in right compounds. " + " Please add it before trying to set the abundance")
        if isinstance(abundance, float):
            self.right_abundance[cmpd] = abundance
        elif isinstance(abundance, int):
            self.right_abundance[cmpd] = float(abundance)
        else:
            raise TypeError("Abundance must be an int or a float")

    def get_right_compound_abundance(self, cmpd):
        """
        Get the abundance of the compound on the right side of the equation.

        :param cmpd: The compound to get the abundance of
        :type cmpd: Compound
        :return: The compounds abundance
        :rtype: float
        """

        if cmpd in self.right_abundance:
            return self.right_abundance[cmpd]
        else:
            raise KeyError(f"In the reaction {self.readable_name} (reaction id: {self.id}), you do not have" +
                           f" {cmpd} on the right hand side of the equation: {self.equation}")

    def number_of_right_compounds(self):
        """
        The number of compounds on the right side of the equation.
        :rtype: int
        """
        return len(self.right_compounds)

    def all_compounds(self):
        """
        Get all the compounds involved in this reaction.
        :return: A set of all the compounds
        :rtype: set
        """
        return self.left_compounds.union(self.right_compounds)

    def number_of_compounds(self):
        """
        Get the total number of compounds involved in this reaction.
        :rtype: int
        """
        return len(self.all_compounds())

    def has(self, cmpd):
        """
        Does this reaction have a compound? Just returns true if the compound is present somewhere in the reaction.
        :param cmpd: The compound to test for
        :type cmpd: Compound
        :rtype: bool
        """
        return cmpd in self.left_compounds or cmpd in self.right_compounds

    def opposite_sides(self, cmpd1, cmpd2):
        """
        Are these two compounds on opposite sides of the reaction?
        :param cmpd1: The first compound
        :type cmpd1: Compound
        :param cmpd2: The second compound
        :type cmpd2: Compound
        :return: Whether the compounds are on opposite sides
        :rtype: bool
        """
        if not self.has(cmpd1):
            raise ValueError(str(cmpd1) + " is not in this reaction")
        if not self.has(cmpd2):
            raise ValueError(str(cmpd2) + " is not in this reaction")
        if cmpd1 in self.left_compounds and cmpd2 in self.right_compounds:
            return True
        if cmpd1 in self.right_compounds and cmpd2 in self.left_compounds:
            return True
        return False

    def set_probability_left_to_right(self, p):
        """
        Set the probability of the reaction running left to right. Note you can also access this as reaction.pLR
        :param p: The probablity
        :type p: float

        """
        if isinstance(p, float):
            self.pLR = p
        elif isinstance(p, int):
            self.pLR = float(p)
        else:
            raise TypeError("The probability must be an int or a float")

    def get_probability_left_to_right(self):
        """
        Get the probability of the reaction running left to right. Note you can also access this as reaction.pLR
        :return: The probablity
        :rtype p: float
        """
        return self.pLR

    def set_probability_right_to_left(self, p):
        """
        Set the probability of the reaction running right to left Note you can also access this as reaction.pRL
        :param p: The probablity
        :type p: float
        """
        if isinstance(p, float):
            self.pRL = p
        elif isinstance(p, int):
            self.pRL = float(p)
        else:
            raise TypeError("The probability must be an int or a float")

    def get_probability_right_to_left(self):
        """
        Get the probability of the reaction running right to left. Note you can also access this as reaction.pRL
        :return: The probablity
        :rtype p: float
        """
        return self.pRL

    def add_enzymes(self, enz):
        """
        Add one or more enzymes that completes this reaction.

        :param enz: A set of enzymes that you want to add
        :type enz: set
        """
        if isinstance(enz, set):
            self.enzymes.update(enz)
        else:
            raise TypeError("You need to supply a set of enzymes")

    def has_enzyme(self, enz):
        """
        Check whether an enzyme is involved in this reaction.
        :param enz: An Enzyme object
        :type enz: Enzyme
        :return: Whether we have this enzyme
        :rtype: bool
        """
        return enz in self.enzymes

    def all_enzymes(self):
        """
        Get all the enzymes involved in this reaction. Returns a set of complex IDs.
        :rtype: set
        """
        return self.enzymes

    def number_of_enzymes(self):
        """
        Gets the number of enzymes involved in this reaction.
        :rtype: int
        """
        return len(self.enzymes)

    def add_pegs(self, pegs):
        """
        Add one or more pegs to this reaction. Pegs must be a set.
        :param pegs: The pegs to add to the reaction
        :type pegs: set
        """
        if isinstance(pegs, set):
            self.pegs.update(pegs)
        else:
            raise TypeError("pegs must be a set")

    def has_peg(self, peg):
        """
        Check whether a peg is involved in this reaction.
        :param peg: The peg to check for
        :type peg: str
        :rtype: bool
        """
        return peg in self.pegs

    def set_deltaG(self, dg):
        """
        Set the value for delta G (Gibbs free energy) for this reaction. Recall -ve deltaG means the reaction is
        favorable.
        :param dg: The delta G of the reaction
        :type dg: float
        """
        if isinstance(dg, float):
            self.deltaG = dg
        elif isinstance(dg, int):
            self.deltaG = float(dg)
        else:
            raise TypeError("The delta G must be an int or a float")
    
    def get_deltaG(self):
        """
        Get the value for delta G (Gibbs free energy) for this reaction.
        :rtype: float
        """
        return self.deltaG

    def check_input_output(self):
        """
        Check whether this reaction is an input or output reaction.
        This is called when we ask is_input_reaction / is_output_reaction and both inp and outp are False
        """

        # do we have external compounds on the left ... then it is an input reaction

        for c in self.left_compounds:
            if c.location == 'e':
                self.inp = True

        for c in self.right_compounds:
            if c.location == 'e':
                self.outp = True

    def toggle_input_reaction(self):
        """
        Set this reaction as an input reaction. This only applies to
        this reaction, so if it is true we set it false, else we set 
        it true
        """
        if self.inp:
            self.inp = False
        else:
            self.inp = True

    def is_input_reaction(self):
        """
        Is this an input reaction?

        :rtype: bool
        """
        if self.inp is False and self.outp is False:
            self.check_input_output()
        return self.inp

    def toggle_output_reaction(self):
        """
        Set this reaction as an output reaction. This only applies to
        this reaction, so if it is true we set it false, else we set 
        it true
        """
        if self.outp:
            self.outp = False
        else:
            self.outp = True

    def is_output_reaction(self):
        """
        Is this an output reaction?
        :rtype: bool
        """
        if self.inp is False and self.outp is False:
            self.check_input_output()
        return self.outp

    def reverse_reaction(self):
        """
        Reverse the reaction - move the left compounds to the right,
        and vice versa. We also switch the abundances and the pLR and
        pRL. 
        
        We also negate the deltaG, since that should be the other way
        around now.
        
        At the moment we don't switch input/output, not sure if we
        need to do that.
        """
        (self.left_compounds, self.right_compounds) = (self.right_compounds, self.left_compounds)
        (self.left_abundance, self.right_abundance) = (self.right_abundance, self.left_abundance)
        (self.inp, self.outp) = (self.outp, self.inp)

        # we only need to reverse two directions
        if self.direction == '>':
            self.direction = '<'
        elif self.direction == '<':
            self.direction = '>'

        # we only need to reverse two gfdirections
        if self.gfdirection == '>':
            self.gfdirection = '<'
        elif self.gfdirection == '<':
            self.gfdirection = '>'

        if self.lower_bound != None and self.upper_bound != None:
            lbtemp = 0 - self.lower_bound
            self.lower_bound = 0 - self.upper_bound
            self.upper_bound = lbtemp

        (self.pLR, self.pRL) = (self.pRL, self.pLR)
        self.deltaG = -self.deltaG

    def add_attribute(self, key, value):
        """
        Add an attribute to this class
        """
        setattr(self, key, value)

    def get_attribute(self, key):
        """
        Retrieve an attribute
        """
        return getattr(self, key)

    def reset_bounds(self):
        """
        reset the bounds of this reaction. If we are using this in gapfilling, we need to reset the bounds
        so we can calculate appropriately.
        :return: None
        """
        self.lower_bound = None
        self.upper_bound = None