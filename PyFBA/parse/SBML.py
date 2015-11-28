import argparse
import os
import sys
from bs4 import BeautifulSoup

import PyFBA

__author__ = 'Rob Edwards'


class SBML:
    """A SBML object representing the model data.

    :ivar model_id: the identifier of the model
    :ivar model_name: the name of the model
    :ivar reactions: a dictionary of Reaction objects with str(obj) as the key
    :ivar compounds: a dictionary of Compound objects with str(obj) as the key
    :ivar compounds_by_id: a dictionary of Compound objects with compound.model_seed_id as the key
    :ivar compartment: a dictionary of compartments in the model

    """

    def __init__(self):
        self.reactions = {}
        self.compounds = {}
        self.compounds_by_id = {}
        self.model_id = ""
        self.model_name = ""
        self.compartment = {}

    def add_compound(self, cpd):
        """
        Add a compound to the model
        :param cpd: The compound to be added as a Compound object
        :type cpd: object.
        :return: None
        :rtype: None
        """
        self.compounds[str(cpd)] = cpd
        if not cpd.abbreviation:
            sys.stderr.write("WARNING: No id for " + cpd.name + "\n")
        else:
            self.compounds_by_id[cpd.abbreviation] = cpd

    def get_all_compounds(self):
        """
        Get the compounds in the model
        :return: A list of the compounds in this model
        :rtype: list of Compound
        """
        return self.compounds.values()

    def get_a_compound(self, cpd):
        """
        Get a single compound with the same str() as the one provided
        :param cpd: The compound to fetch
        :type cpd: object.
        :return: The compound from the model
        :rtype: Compound
        """
        if str(cpd) in self.compounds:
            return self.compounds[str(cpd)]
        else:
            raise ValueError(str(cpd) + " is not present in the model")

    def get_a_compound_by_id(self, cpdid):
        """
        Get a single compound by the id of the compound
        :param cpdid: compound id
        :type cpdid: str
        :return: The compound from the model
        :rtype: Compound
        """

        if cpdid in self.compounds_by_id:
            return self.compounds_by_id[cpdid]
        else:
            raise ValueError(cpdid + " is not present in the model")

    def add_reaction(self, rxn):
        """
        Add a reaction to the model
        :param rxn: The reaction to be added as a metabolism.Reaction object
        :type rxn: object.
        :return: None
        :rtype: None
        """
        self.reactions[str(rxn)] = rxn

    def get_all_reactions(self):
        """
        Get all the reactions
        :return: A set of the reactions in the model
        :rtype: set of Reaction
        """
        return self.reactions

    def get_a_reaction(self, rxn):
        """
        Get a single reaction with the same str() as the one provided
        :param rxn: The reaction to retrieve
        :type rxn: object.
        :return: The reaction object
        :rtype: Reaction

        """

        if (str(rxn)) in self.reactions:
            return self.reactions[str(rxn)]
        else:
            raise ValueError(str(rxn) + " is not present in the model")


def parse_sbml_file(sbml_file, verbose=False):
    """
    Parse an SBML file and return an SBML object.

    :param sbml_file: the SBML file to parse
    :type sbml_file: str
    :param verbose: Whether to create more output
    :type verbose: bool.
    :return: An SBML object
    :rtype: object.
    """

    if not os.path.exists(sbml_file):
        raise IOError("SBML file {} was not found".format(sbml_file))
    soup = BeautifulSoup(open(sbml_file, 'r'), 'xml')
    sbml = SBML()
    sbml.model_name = soup.model['name']
    sbml.model_id = soup.model['id']

    # parse the compartments
    for c in soup.listOfCompartments.find_all('compartment'):
        sbml.compartment[c['id']] = c['name']

    # add the compounds
    for s in soup.listOfSpecies.find_all('species'):
        cpd = PyFBA.metabolism.Compound(s['name'].replace('_c0', '').replace('_e0', ''),
                                        s['compartment'].replace('0', ''))
        cpd.abbreviation = s['id']
        cpd.model_seed_id = s['id'].replace('_c0', '').replace('_e0', '')
        cpd.charge = s['charge']
        if s['boundaryCondition'] == 'false':
            cpd.uptake_secretion = False
        elif s['boundaryCondition'] == 'true':
            cpd.uptake_secretion = True
        else:
            if verbose:
                sys.stderr.write("No boundary rule for {}\n".format(cpd.name))
            cpd.uptake_secretion = False
        sbml.add_compound(cpd)

    # add the reactions
    for r in soup.listOfReactions.find_all('reaction'):
        # I am going to split off the location for the reaction.
        # I don't believe we have the same reaction running in two different locations but maybe in plants, etc?
        if 'biomass' in r['id'].lower():
            rxnid = 'biomass_equation'
        elif '_' not in r['id']:
            if verbose:
                sys.stderr.write("Warning: " + r['id'] + " seems to be a weird id\n")
            rxnid = r['id']
        elif r['id'].startswith('EX_'):
            ex, rxnid, rxnloc = r['id'].split("_")
            rxnid = 'EX_' + rxnid
        else:
            try:
                rxnid, rxnloc = r['id'].split("_")
            except IndexError:
                if verbose:
                    sys.stderr.write("ERROR: Can't unpack " + r['id'] + "\n")
                continue

        rxn = PyFBA.metabolism.Reaction(rxnid)
        if rxn in sbml.get_all_reactions():
            if verbose:
                sys.stderr.write("Already found reaction: " + str(rxn) + " ... not overwriting\n")
            continue
        rxn.description = r['name']
        if rxnid == 'biomass_equation':
            rxn.set_direction('>')
        elif r['reversible'] == 'true':
            rxn.set_direction("=")
        else:
            rxn.set_direction(">")

        # a hash to build the equation from
        equation = {'left': [], 'right': []}
        for rc in r.find_all('listOfReactants'):
            for sp in rc.find_all('speciesReference'):
                cpdname, cpdloc = sp['species'].split("_")
                try:
                    # cpd = sbml.get_a_compound(Compound(cpdname, cpdloc))
                    cpd = sbml.get_a_compound_by_id(sp['species'])
                except ValueError:
                    # the compound is not in the model (but it should be!)
                    cpdnew = PyFBA.metabolism.Compound(cpdname, cpdloc)
                    if verbose:
                        sys.stderr.write("WARNING: {} loc: {}".format(cpdname, cpdloc) +
                                         " is supposed to be in the model but is not. Added\n")
                    sbml.add_compound(cpdnew)
                    cpd = sbml.get_a_compound(PyFBA.metabolism.Compound(cpdname, cpdloc))

                rxn.add_left_compounds({cpd})
                rxn.set_left_compound_abundance(cpd, float(sp['stoichiometry']))
                if cpd.uptake_secretion:
                    rxn.is_uptake_secretion = True
                equation['left'].append(" (" + str(sp['stoichiometry']) + ") " + str(cpd))

        for rc in r.find_all('listOfProducts'):
            for sp in rc.find_all('speciesReference'):
                cpdname, cpdloc = sp['species'].split("_")
                try:
                    # cpd = sbml.get_a_compound(Compound(cpdname, cpdloc))
                    cpd = sbml.get_a_compound_by_id(sp['species'])
                except ValueError:
                    # the compound is not in the model (but it should be!)
                    cpdnew = PyFBA.metabolism.Compound(cpdname, cpdloc)
                    if verbose:
                        sys.stderr.write("WARNING: {} is supposed to be in the model but is not.\n".format(cpdnew))
                    sbml.add_compound(cpdnew)
                    cpd = sbml.get_a_compound(PyFBA.metabolism.Compound(cpdname, cpdloc))

                rxn.add_right_compounds({cpd})
                rxn.set_right_compound_abundance(cpd, float(sp['stoichiometry']))
                if cpd.uptake_secretion:
                    rxn.is_uptake_secretion = True
                equation['right'].append(" (" + str(sp['stoichiometry']) + ") " + str(cpd))

        rxn.equation = " + ".join(equation['left']) + " " + rxn.direction + " " + " + ".join(equation['right'])

        for params in r.find_all('listOfParameters'):
            for p in params.find_all('parameter'):
                if p['id'].lower() == 'lower_bound':
                    rxn.lower_bound = float(p['value'])
                if p['id'].lower() == 'upper_bound':
                    rxn.upper_bound = float(p['value'])

        sbml.add_reaction(rxn)

    return sbml

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse an SBML file")
    parser.add_argument('-s', help='SBML file', required=True)
    args = parser.parse_args()

    sbml_data = parse_sbml_file(args.s, True)
    print("There are " + str(len(sbml_data.get_all_reactions())) + " reactions and " +
          str(len(sbml_data.get_all_compounds())) + " compounds present in the model\n")
