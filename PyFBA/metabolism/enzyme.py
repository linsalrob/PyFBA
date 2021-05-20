"""
The enzyme is responsible for running a reaction. It also connects our
reactions to a genome. We need to know which enzymes we have in the
genome we are investigating

"""
from . import Reaction


class Enzyme:
    """
    The enzyme class has a few components:
      * The subunit(s) that make up the enzyme
      * The genes that encode those subunit(s)
      * The reactions that this enzyme is connected to.

    :ivar name: the name of the enzyme object
    :type name: str
    :ivar roles: the set of roles associated with the enzyme object
    :type roles: set
    :ivar pegs: a dict of pegs associated with the enzyme object and their associated roles
    :type pegs: dict
    :ivar roles_w_pegs: a dict of roles associated with the enzyme and their pegs
    :type roles_w_pegs: dict
    :ivar reactions: a set of reaction IDs that this enzyme connects to
    :type reactions: set
    :ivar ec_number: one or more EC numbers associated with this Enzyme. We only store the numeric part (not "EC: ")
    :type ec_number: set
    """

    def __init__(self, name):
        """
        Instantiate the enzyme

        :param name: the name of the enzyme
        :type name: str
        """

        self.name = name  # whatever name we give to this thing!
        self.roles = set()  # Roles (text strings)
        self.pegs = {}  # a hash that connects Roles to PEGs
        self.roles_w_pegs = {}  # which roles have pegs
        self.reactions = set()  # RIDs that the enzyme connects to
        self.ec_number = set()  # one or more EC numbers associated with this Enzyme.

    def __eq__(self, other):
        """
        Is this enzyme the same as another one?

        :param other: The other enzyme
        :type other: Enzyme
        :return: Whether the two enzymes are the same
        :rtype: bool
        """
        if isinstance(other, Enzyme):
            return (self.name, self.roles) == (other.name, other.roles)
        else:
            return NotImplemented

    def __ne__(self, other):
        """
        Are these not equal?

        :param other: The other enzyme
        :type other: Enzyme
        :return: Whether the two enzymes are not equal
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
        return hash(self.name)

    def __str__(self):
        """
        The string representation of the enzyme

        :rtype: str
        """
        return "ENZYME: " + self.name + " (roles: " + "; ".join([x for x in self.roles]) + ")"

    def add_roles(self, roles):
        """
        Add roles to this enzyme or complex

        :param roles: A set of functional roles that encode the enzyme
        :type roles: set
        """
        if not isinstance(roles, set):
            raise TypeError("Roles must be a set")
        self.roles.update(roles)
    
    def has_role(self, role):
        """
        Does this enzyme have this role?

        :param role: A functional role
        :type role: str
        :returns: A boolean
        :rtype: bool
        """
        return role in self.roles

    def number_of_roles(self):
        """
        How many roles does this enzyme have?

        :rtype: int
        """
        return len(self.roles)

    def add_pegs(self, pegs):
        """
        Add a hash of pegs and roles. Keys must be pegs, values must be roles.

        Will throw a KeyError if the Role is not present

        :param pegs: A hash of pegs and roles that encode the enzyme (e.g. from the assigned functions file)
        :type pegs: dict
        :raises: KeyError
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
        """
        Add a single peg and the role that it connects to.

        :param peg: The peg id
        :type peg: str
        :param role: The role it connects to
        :type role: str
        :raises: KeyError
        """
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
        """
        The number of pegs associated with this enzyme.

        :rtype: int
        """
        return len(self.pegs)

    def number_of_roles_with_pegs(self):
        """
        How many of our roles have pegs associated with them?

        :rtype: int
        """
        return len(self.roles_w_pegs)

    def has_peg_for_role(self, role):
        """
        Do we have at least one peg for this role?

        :param role: The role we are looking for
        :type role: str
        :return: If a peg is present
        :rtype: bool
        """
        return role in self.roles_w_pegs

    def add_reaction(self, reaction):
        """
        Add a reaction that this enzyme is inolved in.

        :param reaction: The reaction object that this is involved with
        :type reaction: Reaction
        """
        if not isinstance(reaction, str):
            raise TypeError("reaction must be a string not a " + str(type(reaction)))
        self.reactions.add(reaction)

    def number_of_reactions(self):
        """
        The number of reactions that this enzyme is involved in

        :rtype: int
        """
        return len(self.reactions)

    def add_ec(self, ecnumber):
        """
        Add an EC number to the Enzyme complex. We just store the 1.2.3.4 or 1.2.-.- part, not the EC part.

        :param ecnumber: The EC number
        :type ecnumber: str
        """

        self.ec_number.add(ecnumber)

    def probability(self):
        """
        The probability that this reaction occurs in the cell.
        Currently the number of pegs/number of roles. Thus if most of
        the pegs are present then the enzyme is likely to function

        :returns: the probability that this reaction is complete
        :rtype: float
        """
        
        # Initially we had this, but note that a peg can have two roles
        # (joined together with " / " or " @ ", and so we can't just use
        # this simple calculation. We need to know thenumber of pegroles
        # / number of roles - roles with pegs!
        return 1.0 * self.number_of_roles_with_pegs() / self.number_of_roles()
