import argparse
import copy
import os
import sys
from bs4 import BeautifulSoup

import PyFBA
from PyFBA import log_and_message


class SBML:
    """A SBML object representing the model data.

    :ivar model_id: the identifier of the model
    :ivar model_name: the name of the model
    :ivar reactions: a dictionary of Reaction objects with id as the key
    :ivar compounds: a dictionary of Compound objects with id as the key
    :ivar compounds_by_name: a dictionary of Compound objects with compound.model_seed_id as the key
    :ivar compartment: a dictionary of compartments in the model

    """

    def __init__(self):
        self.reactions = {}
        self.compounds = set()
        self.compounds_by_id = {}
        self.compounds_by_name = {}
        self.compounds_by_id_loc = {}
        self.compounds_by_name_loc = {}
        self.compounds_by_id_loc = {}
        self.model_id = ""
        self.model_name = ""
        self.compartment = {}

    def add_compound(self, cpd):
        """
        Add a compound to the model
        :param cpd: The compound to be added as a Compound object
        :type cpd: PyFBA.metabolism.CompoundWithLocation
        :return: None
        :rtype: None
        """

        self.compounds.add(cpd)
        if cpd.id not in self.compounds_by_id:
            self.compounds_by_id[cpd.id] = set()
        if cpd.name not in self.compounds_by_name:
            self.compounds_by_name[cpd.name] = set()
        self.compounds_by_id[cpd.id].add(cpd)
        self.compounds_by_name[cpd.name].add(cpd)

    def get_all_compounds(self):
        """
        Get the compounds in the model
        :return: A set of the compounds in this model
        :rtype: set[PyFBA.metabolism.CompoundWithLocation]
        """
        return self.compounds

    def get_a_compound_by_id(self, cpdid) -> PyFBA.metabolism.CompoundWithLocation:
        """
        Get a single compound by the id of the compound
        :param cpdid: compound id
        :type cpdid: str
        :return: The compound from the model
        :rtype: PyFBA.metabolism.CompoundWithLocation
        """

        if cpdid in self.compounds_by_id:
            return next(iter(self.compounds_by_id[cpdid]))
        # raise ValueError(cpdid + " is not present in the model")
        return None

    def get_a_compound_by_name(self, name) -> PyFBA.metabolism.CompoundWithLocation:
        """
        Get a single compound by its name
        :param name: the name of the compount
        :type name: str
        :return: the compound object
        :rtype: PyFBA.metabolism.CompoundWithLocation
        """

        if name in self.compounds_by_name:
            return next(iter(self.compounds_by_name[name]))
        # raise ValueError(name + " is not present in the model")
        return None

    def get_a_compound_by_id_and_loc(self, cpdid, loc) -> PyFBA.metabolism.CompoundWithLocation:
        """
        Get a single compound by the id of the compound
        :param cpdid: compound id
        :type cpdid: str
        :param loc: the location
        :type loc: str
        :return: The compound from the model
        :rtype: PyFBA.metabolism.CompoundWithLocation
        """

        if cpdid in self.compounds_by_id:
            for c in self.compounds_by_id[cpdid]:
                if c.location == loc:
                    return c
        # raise ValueError(f"{cpdid} with location {loc} is not present in the model")
        return None

    def get_a_compound_by_name_and_loc(self, name, loc) -> PyFBA.metabolism.CompoundWithLocation:
        """
        Get a single compound by its name
        :param name: the name of the compount
        :type name: str
        :param loc: the location
        :type loc: str
        :return: the compound object
        :rtype: PyFBA.metabolism.CompoundWithLocation
        """

        if name in self.compounds_by_name:
            for c in self.compounds_by_name[name]:
                if c.location == loc:
                    return c
        # raise ValueError(f"{name} with location {loc} is not present in the model")
        return None

    def add_reaction(self, rxn):
        """
        Add a reaction to the model
        :param rxn: The reaction to be added as a metabolism.Reaction object
        :type rxn: object.
        :return: None
        :rtype: None
        """
        self.reactions[rxn.id] = rxn

    def get_all_reactions(self):
        """
        Get all the reactions
        :return: A set of the reactions in the model
        :rtype: set of Reaction
        """
        return set(self.reactions.values())

    def get_a_reaction(self, rid):
        """
        Get a single reaction with the same id as the one provided
        :param rid: The reaction to retrieve
        :type rid: object.
        :return: The reaction object
        :rtype: Reaction
        """

        if rid in self.reactions:
            return self.reactions[rid]
        else:
            raise ValueError(f"{rid} is not present in the model")

    def correct_media(self, media):
        """
        Correct media components for anything added from the SBML file. Note we have a separate method
        for this that will also work with a set of compounds
        :param media: a set of media compounds
        :type media: set(PyFBA.metabolism.CompoundWithLocation)
        :return: a set of media compounds with corrected names
        :rtype: set(PyFBA.metabolism.CompoundWithLocation)
        """

        new_media = set()
        warned_compounds = False
        for m in media:
            comp = self.get_a_compound_by_name_and_loc(m.name, 'e')
            if comp:
                new_media.add(comp)
                continue

            comp = self.get_a_compound_by_name(m.name)
            if comp:
                media_component = PyFBA.metabolism.CompoundWithLocation.from_compound(comp, 'e')
                new_media.add(media_component)
                continue

            testname = m.name.replace('-', '_')
            comp = self.get_a_compound_by_name_and_loc(testname, 'e')
            if comp:
                new_media.add(self.get_a_compound_by_name_and_loc(testname, 'e'))
                continue

            comp = self.get_a_compound_by_name(testname)
            if comp:
                media_component = PyFBA.metabolism.CompoundWithLocation.from_compound(comp, 'e')
                new_media.add(media_component)
                continue

            testname = m.name.replace('+', '')
            comp = self.get_a_compound_by_name_and_loc(testname, 'e')
            if comp:
                new_media.add(self.get_a_compound_by_name_and_loc(testname, 'e'))
                continue

            comp = self.get_a_compound_by_name(testname)
            if comp:
                media_component = PyFBA.metabolism.CompoundWithLocation.from_compound(comp, 'e')
                new_media.add(media_component)
                continue

            log_and_message(f"Checking media compounds: Our compounds do not include  {m.name}", stderr=True)
            warned_compounds = True
            new_media.add(m)

        if warned_compounds:
            log_and_message("""
Please note: The warnings about media not being found in compounds are not fatal.
It just means that we did not find that compound anywhere in the reactions, and so it is unlikely to be
needed or used. We typically see a few of these in rich media. 
            """, stderr=True)
        return new_media


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
        if 'name' in c:
            sbml.compartment[c['id']] = c['name']
        else:
            sbml.compartment[c['id']] = c['id']

    # add the compounds
    for s in soup.listOfSpecies.find_all('species'):
        cpdid = s['id'].replace('_c0', '').replace('_e0', '')
        cpdname = s['name'].replace('_c0', '').replace('_e0', '')
        cpdloc = s['compartment'].replace('0', '')
        if cpdid.startswith('M_'):
            cpdid = cpdid.replace('M_', "")
        if '_b' in cpdid:
            cpdloc = 'b'
            cpdid = cpdid.replace('_b', '')
            cpdname = cpdname.replace('_b', '')
        cpd = PyFBA.metabolism.CompoundWithLocation(cpdid, cpdname, cpdloc)
        cpd.abbreviation = s['id']
        cpd.model_seed_id = cpdid
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
        if 'biomass' in r['id'].lower():
            rxnid = 'biomass_equation'
        elif '_' not in r['id']:
            if verbose:
                sys.stderr.write("Warning: " + r['id'] + " seems to be a weird id\n")
            rxnid = r['id']
        elif r['id'].startswith('EX_'):
            ex, rxnid, rxnloc = r['id'].split("_")
            rxnid = 'EX_' + rxnid
        elif r['id'].startswith('R_'):
            ex, rxnid, rxnloc = r['id'].split("_")
        else:
            try:
                rxnid, rxnloc = r['id'].split("_")
            except (ValueError, IndexError):
                if verbose:
                    sys.stderr.write("ERROR: Can't unpack " + r['id'] + "\n")
                continue

        rxn = PyFBA.metabolism.Reaction(rxnid)

        if rxn in sbml.get_all_reactions():
            log_and_message("Already found reaction: {rxn} ... not overwriting", stderr=verbose)
            continue
        rxn.readable_name = r['name']
        if rxnid == 'biomass_equation':
            rxn.set_direction('>')
        elif r['reversible'] == 'true':
            rxn.set_direction("=")
        else:
            rxn.set_direction(">")

        # a hash to build the equation from
        equation = {'left': [], 'right': []}
        count = 0

        # here we find the reactants and products from the SBML file and
        # add them to the left and right equation arrays appropriately.

        for rp in ['listOfReactants', 'listOfProducts']:
            for rc in r.find_all(rp):
                for sp in rc.find_all('speciesReference'):
                    if sp['species'].startswith('M_'):
                        m, cpdid, cpdloc = sp['species'].split("_")
                    else:
                        cpdid, cpdloc = sp['species'].split("_")
                    cpdloc = cpdloc.replace('0', '')

                    try:
                        cpd = copy.deepcopy(sbml.get_a_compound_by_id_and_loc(cpdid, cpdloc))
                    except ValueError:
                        # the compound is not in the model (but it should be!)
                        count += 1
                        cpd = PyFBA.metabolism.CompoundWithLocation(f"smbl{count}", cpdid, cpdloc)
                        log_and_message(
                            f"WARNING: {cpdid} loc: {cpdloc} is supposed to be in the model but is not. Added\n",
                            c="RED",
                            stderr=verbose
                        )
                        sbml.add_compound(cpd)

                    if isinstance(cpd, PyFBA.metabolism.CompoundWithLocation):
                        cpd.location = cpdloc
                    elif isinstance(cpd, PyFBA.metabolism.Compound):
                        cpd = PyFBA.metabolism.CompoundWithLocation.from_compound(cpd, cpdloc)
                    else:
                        log_and_message(f"ERROR: {cpd} is neither Compound nor CompoundWithLocation", stderr=verbose)

                    if cpd.uptake_secretion:
                        rxn.is_uptake_secretion = True

                    if 'listOfReactants' == rp:
                        rxn.add_left_compounds({cpd})
                        rxn.set_left_compound_abundance(cpd, float(sp['stoichiometry']))
                        equation['left'].append(" (" + str(sp['stoichiometry']) + ") " + str(cpd))
                    else:
                        rxn.add_right_compounds({cpd})
                        rxn.set_right_compound_abundance(cpd, float(sp['stoichiometry']))
                        equation['right'].append(" (" + str(sp['stoichiometry']) + ") " + str(cpd))

        rxn.equation = " + ".join(equation['left']) + " " + rxn.direction + " " + " + ".join(equation['right'])

        for params in r.find_all('listOfParameters'):
            for p in params.find_all('parameter'):
                if p['id'].lower() == 'lower_bound':
                    rxn.lower_bound = float(p['value'])
                if p['id'].lower() == 'upper_bound':
                    rxn.upper_bound = float(p['value'])

        sbml.add_reaction(rxn)

    log_and_message(f"Parsing the model {sbml.model_name} (id {sbml.model_id}) is complete.")
    log_and_message(f"Parsing the SBML file: We found {len(sbml.compounds)} compounds")
    log_and_message(f"Parsing the SBML file: We found {len(sbml.reactions)} reactions")

    return sbml



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse an SBML file")
    parser.add_argument('-s', help='SBML file', required=True)
    args = parser.parse_args()

    sbml_data = parse_sbml_file(args.s, True)
    print("There are " + str(len(sbml_data.get_all_reactions())) + " reactions and " +
          str(len(sbml_data.get_all_compounds())) + " compounds present in the model\n")
