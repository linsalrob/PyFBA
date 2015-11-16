"""
The enzyme is responsible for running a reaction. It also connects our
reactions to a genome. We need to know which enzymes we have in the
genome we are investigating

"""


class Enzyme():
    """
    The enzyme class has a few components
    
    The subunit(s) that make up the enzyme

    The genes that encode those subunit(s)

    The reactions that this enzyme is connected to.

    """


    def __init__(self, name):
        self.name = name # whatever name we give to this thing!
        self.roles = set() # Roles (text strings)
        self.pegs = {} # a hash that connects Roles to PEGs
        self.roles_w_pegs = {} # which roles have pegs
        self.reactions = set() # RIDs that the enzymeconnects to

    def __eq__(self, other):
        """ is this enzyme the same as another one? """
        if isinstance(other, Enzyme):
            return ((self.name, self.roles) == (other.name, other.roles))
        else:
            return NotImplemented

    def __ne__(self, other):
        """Are these not equal?"""
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __hash__(self):
        """The hash function is based on the name of the compound"""
        return hash(self.name)

    def __str__(self):
        """The to string function"""
        return "ENZYME: " + self.name + " (roles: " + "; ".join([x for x in self.roles]) + ")"

    def add_roles(self, roles):
        """Add roles to this enzyme or complex"""
        if not isinstance(roles, set):
            raise TypeError("Roles must be a set")
        self.roles.update(roles)
    
    def has_role(self, role):
        """Does this enzyme have this role"""
        return role in self.roles

    def number_of_roles(self):
        """How many roles does this enzyme have?"""
        return len(self.roles)

    def add_pegs(self, pegs):
        """Add a hash of pegs and roles. Keys must be pegs, values must
        be roles.

        Will throw a KeyError if the Role is not present
        """
        if not isinstance(pegs, dict):
            raise TypeError("pegs must be a hash to add more than one")
        for p in pegs:
            if pegs[p] in self.roles:
                self.pegs[p] = pegs[p]
                if pegs[p] not in self.roles_w_pegs:
                    self.roles_w_pegs[pegs[p]] = []
                self.roles_w_pegs[pegs[p]].append(p)
            else:
                raise KeyError("Role " + pegs[p] + " not found")

    def add_a_peg(self, peg, role):
        """Add a single peg and the role that it connects to"""
        if not isinstance(peg, str):
            raise TypeError("peg must be a string. Did you mean to use add_pegs?")
        if role in self.roles:
            self.pegs[peg] = role
            if role not in self.roles_w_pegs:
                self.roles_w_pegs[role] = []
            self.roles_w_pegs[role].append(peg)
        else:
            raise KeyError("Role " + role + " not found")

    def number_of_pegs(self):
        """The number of pegs assocaited with this enzyme"""
        return len(self.pegs)

    def number_of_roles_with_pegs(self):
        """How many of our roles have pegs associated with them?"""
        return len(self.roles_w_pegs)

    def has_peg_for_role(self, role):
        """Do we have at least one peg for this role?"""
        return role in self.roles_w_pegs

    def add_reaction(self, reaction):
        """Add a reaction that this enzyme is inolved in"""
        if not isinstance(reaction, str):
            raise TypeError("reaction must be a string")
        self.reactions.add(reaction)

    def number_of_reactions(self):
        """The number of reactions that this enzyme is involved in"""
        return len(self.reactions)


    def probability(self):
        """The probability that this reaction occurs in the cell.
        Currently the number of pegs/number of roles. Thus if most of
        the pegs are present then the enzyme is likely to function
        """
        
        # Initially we had this, but note that a peg can have two roles
        # (joined together with " / " or " @ ", and so we can't just use
        # this simple calculation. We need to know thenumber of pegroles
        # / number of roles - roles with pegs!
        #return 1.0 * self.number_of_pegs()/self.number_of_roles()
        return 1.0 * self.number_of_roles_with_pegs() / self.number_of_roles()


