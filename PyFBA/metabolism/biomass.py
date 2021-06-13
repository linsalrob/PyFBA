import sys
import PyFBA
from .compound import CompoundWithLocation
from .reaction import Reaction

def standard_eqn():
    """The standard biomass_equation equation is derived from the SBML file
    (actually for C. sedlakii, but I suspect it is the same for others)
    and includes amino acids, nucleotides.

    :returns: the reactants and products
    :rtype: dict, dict
    """

    reactants = {
        "ATP": 40.110, "L-Valine": 0.353, "L-Alanine": 0.429, "L-Methionine": 0.128, "L-Phenylalanine": 0.155,
        "S-Adenosyl-L-methionine": 0.008, "CoA": 0.008, "CTP": 0.084, "TTP": 0.015, "dCTP": 0.015,
        "L-Isoleucine": 0.242, "CO2": 0.008, "L-Aspartate": 0.201, "L-Glutamate": 0.219, "L-Glutamine": 0.219,
        "L-Leucine": 0.376, "L-Serine": 0.180, "L-Tryptophan": 0.047, "L-Proline": 0.185, "GTP": 0.135, "dGTP": 0.015,
        "L-Threonine": 0.211, "NAD": 0.008, "dATP": 0.015, "L-Tyrosine": 0.115, "L-Asparagine": 0.201, "NADP": 0.008,
        "L-Cysteine": 0.076, "L-Histidine": 0.079, "L-Lysine": 0.286, "UTP": 0.091, "L-Arginine": 0.247,
        "Glycine": 0.511
    }

    products = {
        "ADP": 40.000, "Phosphate": 39.992, "apo-ACP": 0.008, "Biomass": 1.000, "H": 40.000, "PPi ": 0.406
    }

    return reactants, products


def kbase():
    """The kbase biomass_equation reaction.

    :returns: the reactants and products
    :rtype: dict, dict
    """
    reactants = {'GTP': 0.135406821203723, 'L-Aspartate': 0.200830806928348, 'L-Methionine': 0.127801422590767,
                 'CTP': 0.0841036156544863, 'NAD': 0.00309646685192537, 'Fe3+': 0.00309646685192537,
                 'K+': 0.00309646685192537, 'DNA replication': 1.0, 'Sulfate': 0.00309646685192537,
                 'dATP': 0.0145080770930701, 'Ubiquinone-8': 0.00309646685192537, 'ACP': 0.00309646685192537,
                 'Mn2+': 0.00309646685192537, 'ATP': 40.1101757365074, 'GSH': 0.00309646685192537,
                 'Zn2+': 0.00309646685192537, 'L-Threonine': 0.211072732780569, 'Heme': 0.00309646685192537,
                 'L-Glutamate': 0.219088153012743, 'Phosphatidylglycerol dioctadecanoyl': 0.0106480421341882,
                 'Protein biosynthesis': 1.0, 'Peptidoglycan polymer (n subunits)': 0.0250105977108944,
                 'L-Lysine': 0.285438020490179, 'phosphatidylethanolamine dioctadecanoyl': 0.0106480421341882,
                 'Siroheme': 0.00309646685192537, 'L-Asparagine': 0.200830806928348, 'L-Valine': 0.352233189091625,
                 'TTP': 0.0145080770930701, 'Putrescine': 0.00309646685192537, 'Co2+': 0.00309646685192537,
                 'Spermidine': 0.00309646685192537, 'Ca2+': 0.00309646685192537, 'L-Phenylalanine': 0.154519490031345,
                 'L-Leucine': 0.375388847540127, 'L-Cysteine': 0.0761464922056484,
                 '2-Demethylmenaquinone 8': 0.00309646685192537, 'L-Serine': 0.179456352975885,
                 'dCTP': 0.017531703978307, '10-Formyltetrahydrofolate': 0.00309646685192537,
                 'L-Glutamine': 0.219088153012743, 'FAD': 0.00309646685192537,
                 'Bactoprenyl diphosphate': 0.0250105977108944, 'Cl-': 0.00309646685192537,
                 'Mg': 0.00309646685192537, 'Anteisoheptadecanoylcardiolipin (B. subtilis)': 0.0106480421341882,
                 'L-Tyrosine': 0.120676604606612, 'Dianteisoheptadecanoylphosphatidylglycerol': 0.0106480421341882,
                 'S-Adenosyl-L-methionine': 0.00309646685192537, 'L-Histidine': 0.0792636000737159,
                 'L-Arginine': 0.246696822701341, 'Menaquinone 8': 0.00309646685192537,
                 '5-Methyltetrahydrofolate': 0.00309646685192537, 'UTP': 0.0908319049068452,
                 'Isoheptadecanoylcardiolipin (B. subtilis)': 0.0106480421341882, 'TPP': 0.00309646685192537,
                 'NADP': 0.00309646685192537, 'L-Isoleucine': 0.241798510337235,
                 'Diisoheptadecanoylphosphatidylethanolamine': 0.0106480421341882, 'L-Alanine': 0.427934380173264,
                 'Stearoylcardiolipin (B. subtilis)': 0.0106480421341882,
                 'core oligosaccharide lipid A': 0.0250105977108944, 'Glycine': 0.509869786991038,
                 'Tetrahydrofolate': 0.00309646685192537, 'Calomide': 0.00309646685192537, 'Fe2+': 0.00309646685192537,
                 'dGTP': 0.017531703978307, 'Riboflavin': 0.00309646685192537, 'CoA': 0.00309646685192537,
                 'Dianteisoheptadecanoylphosphatidylethanolamine': 0.0106480421341882, 'H2O': 35.5403092430435,
                 'Pyridoxal phosphate': 0.00309646685192537, 'L-Tryptophan': 0.0472019191450218,
                 'RNA transcription': 1.0, 'L-Proline': 0.184354665339991,
                 'Diisoheptadecanoylphosphatidylglycerol': 0.0106480421341882, 'Cu2+': 0.00309646685192537}
    products = {'apo-ACP': 0.00309646685192537, 'Peptidoglycan polymer (n-1 subunits)': 0.0250105977108944,
                'Dimethylbenzimidazole': 0.00309646685192537, 'PPi': 0.484597640415186, 'Biomass': 1.0, 'H+': 40.0,
                'Phosphate': 39.9969035331481, 'ADP': 40.0, 'Cobinamide': 0.00309646685192537}

    return reactants, products


def kbase_simple():
    """A slightly simplified version of the kbase biomass_equation reaction

    :returns: the reactants and products
    :rtype: dict, dict
    """
    reactants = {'GTP': 0.135406821203723, 'L-Aspartate': 0.200830806928348, 'L-Methionine': 0.127801422590767,
                 'CTP': 0.0841036156544863, 'NAD': 0.00309646685192537, 'Fe3+': 0.00309646685192537,
                 'K+': 0.00309646685192537,
                 'Sulfate': 0.00309646685192537, 'dATP': 0.0145080770930701, 'Ubiquinone-8': 0.00309646685192537,
                 'ACP': 0.00309646685192537, 'Mn2+': 0.00309646685192537, 'ATP': 40.1101757365074,
                 'GSH': 0.00309646685192537,
                 'Zn2+': 0.00309646685192537, 'L-Threonine': 0.211072732780569, 'Heme': 0.00309646685192537,
                 'L-Glutamate': 0.219088153012743, 'Phosphatidylglycerol dioctadecanoyl': 0.0106480421341882,
                 'Peptidoglycan polymer (n subunits)': 0.0250105977108944, 'L-Lysine': 0.285438020490179,
                 'phosphatidylethanolamine dioctadecanoyl': 0.0106480421341882, 'Siroheme': 0.00309646685192537,
                 'L-Asparagine': 0.200830806928348, 'L-Valine': 0.352233189091625, 'TTP': 0.0145080770930701,
                 'Putrescine': 0.00309646685192537, 'Co2+': 0.00309646685192537, 'Spermidine': 0.00309646685192537,
                 'Ca2+': 0.00309646685192537, 'L-Phenylalanine': 0.154519490031345, 'L-Leucine': 0.375388847540127,
                 'L-Cysteine': 0.0761464922056484, '2-Demethylmenaquinone 8': 0.00309646685192537,
                 'L-Serine': 0.179456352975885,
                 'dCTP': 0.017531703978307, '10-Formyltetrahydrofolate': 0.00309646685192537,
                 'L-Glutamine': 0.219088153012743,
                 'FAD': 0.00309646685192537, 'Bactoprenyl diphosphate': 0.0250105977108944, 'Cl-': 0.00309646685192537,
                 'Mg': 0.00309646685192537, 'Anteisoheptadecanoylcardiolipin (B. subtilis)': 0.0106480421341882,
                 'L-Tyrosine': 0.120676604606612, 'Dianteisoheptadecanoylphosphatidylglycerol': 0.0106480421341882,
                 'S-Adenosyl-L-methionine': 0.00309646685192537, 'L-Histidine': 0.0792636000737159,
                 'L-Arginine': 0.246696822701341,
                 'Menaquinone 8': 0.00309646685192537, '5-Methyltetrahydrofolate': 0.00309646685192537,
                 'UTP': 0.0908319049068452,
                 'Isoheptadecanoylcardiolipin (B. subtilis)': 0.0106480421341882, 'TPP': 0.00309646685192537,
                 'NADP': 0.00309646685192537,
                 'L-Isoleucine': 0.241798510337235, 'Diisoheptadecanoylphosphatidylethanolamine': 0.0106480421341882,
                 'L-Alanine': 0.427934380173264, 'Stearoylcardiolipin (B. subtilis)': 0.0106480421341882,
                 'core oligosaccharide lipid A': 0.0250105977108944, 'Glycine': 0.509869786991038,
                 'Tetrahydrofolate': 0.00309646685192537,
                 'Calomide': 0.00309646685192537, 'Fe2+': 0.00309646685192537, 'dGTP': 0.017531703978307,
                 'Riboflavin': 0.00309646685192537,
                 'CoA': 0.00309646685192537, 'Dianteisoheptadecanoylphosphatidylethanolamine': 0.0106480421341882,
                 'H2O': 35.5403092430435, 'Pyridoxal phosphate': 0.00309646685192537,
                 'L-Tryptophan': 0.0472019191450218,
                 'L-Proline': 0.184354665339991, 'Diisoheptadecanoylphosphatidylglycerol': 0.0106480421341882,
                 'Cu2+': 0.00309646685192537}
    products = {'apo-ACP': 0.00309646685192537, 'Peptidoglycan polymer (n-1 subunits)': 0.0250105977108944,
                'Dimethylbenzimidazole': 0.00309646685192537,
                'PPi': 0.484597640415186, 'Biomass': 1.0, 'H+': 40.0, 'Phosphate': 39.9969035331481, 'ADP': 40.0,
                'Cobinamide': 0.00309646685192537}

    return reactants, products


def gram_negative():
    """
    The model from a gap filled Gram_negative bacteria.

    :returns: the reactants and products
    :rtype: dict, dict
    """

    reactants = {
        'K+': 0.00778132482043, 'Mg': 0.00778132482043, 'H2O': 35.5386858538, 'Sulfate': 0.00778132482043,
        'ATP': 40.1101757365, 'RNA transcription': 1, 'L-Phenylalanine': 0.154807600875,
        'L-Serine': 0.179790960094, 'Protein biosynthesis': 1, 'L-Cysteine': 0.0762884719009,
        'NAD': 0.00778132482043, 'L-Arginine': 0.247156803702, 'S-Adenosyl-L-methionine': 0.00778132482043,
        'CoA': 0.00778132482043, 'Fe2+': 0.00778132482043, 'L-Asparagine': 0.201205267996,
        'Peptidoglycan polymer (n subunits)': 0.0609084652443, 'ACP': 0.00778132482043, 'CTP': 0.0841036156545,
        'dATP': 0.0146849834202, 'Ca2+': 0.00778132482043, 'Pyridoxal phosphate': 0.00778132482043,
        'L-Histidine': 0.0794113918032, 'dGTP': 0.0146849834202, 'L-Leucine': 0.376088782529,
        'L-Valine': 0.352889948968, 'L-Alanine': 0.428732289454, 'Cl-': 0.00778132482043, 'DNA replication': 1,
        'UTP': 0.0908319049068, 'Cu2+': 0.00778132482043, 'L-Isoleucine': 0.242249358141,
        'L-Glutamate': 0.219496655995, 'L-Tryptophan': 0.0472899299502, 'Zn2+': 0.00778132482043,
        'Riboflavin': 0.00778132482043, 'TTP': 0.0146849834202, 'Glycine': 0.510820469745, 'Co2+': 0.00778132482043,
        'L-Lysine': 0.285970236775, 'L-Methionine': 0.128039715997, 'FAD': 0.00778132482043,
        'Bactoprenyl diphosphate': 0.0609084652443, 'dCTP': 0.0146849834202, 'L-Tyrosine': 0.115101904973,
        'Mn2+': 0.00778132482043, 'L-Aspartate': 0.201205267996, 'Fe3+': 0.00778132482043,
        'L-Glutamine': 0.219496655995, 'GTP': 0.135406821204, 'NADP': 0.00778132482043,
        'L-Threonine': 0.211466290532, 'L-Proline': 0.184698405655,
    }
    products = {
        'ADP': 40, 'Biomass': 1, 'apo-ACP': 0.00778132482043, 'Phosphate': 39.9922186752, 'PPi': 0.405833094852,
        'Peptidoglycan polymer (n-1 subunits)': 0.0609084652443, 'H+': 40,
    }

    return reactants, products


def biomass_equation(biomass_type='standard', cpds=None, verbose=False):
    """Get the biomass_equation equation for a specific type of biomass_equation equation.

    biomass_type can be one of:
        standard:       the standard biomass_equation equation we were using for the JSON models initially
        kbase:          the revised biomass_equation equation that was included in the kbase models
        kbase_simple:   a simplified version of the kbase biomass_equation equation
        gram_negative:  a Gram negative biomass_equation equation

    compounds is the compounds set Set[PyFBA.metabolism.Compound] that we use to map to the compounds we
    extract from the biomass equations. If none, we will use the modelseed compounds
    :param biomass_type: The type of biomass_equation equation to get
    :type biomass_type: str
    :param cpds: The metabolism compound set that we use to map to the compounds we find
    :type cpds: Set[PyFBA.metabolism.Compound]
    :return: The biomass_equation equation as a Reaction object
    :rtype: Reaction
    """

    modelseed = None

    if isinstance(cpds, set):
        modelseed = PyFBA.model_seed.ModelData(compounds=cpds)
    elif not cpds:
        modelseed = PyFBA.parse.parse_model_seed_data()
    else:
        PyFBA.log_and_message(f"biomass.py can't parse compounds {cpds}", stderr=True)
        return None

    if biomass_type == 'standard':
        reactants, products = standard_eqn()
    elif biomass_type == 'kbase':
        reactants, products = kbase()
    elif biomass_type == 'kbase_simple':
        reactants, products = kbase_simple()
    elif biomass_type == 'gram_negative' or biomass_type == 'gramnegative':
        reactants, products = gram_negative()
    else:
        sys.exit("ERROR: Do not understand what " + biomass_type + " is for a biomass_equation equation\n")

    r = Reaction('biomass_equation', 'biomass_equation')
    for i,c in enumerate(reactants):
        cpdname = modelseed.get_compound_by_name(c)
        if cpdname:
            cpd = CompoundWithLocation.from_compound(cpdname, 'c')
        else:
            PyFBA.log_and_message(f"biomass.py: No compound found for {c} in our compounds dataset", stderr=verbose)
            cpd = CompoundWithLocation(f"rctn{i}", c, 'c')
        r.add_left_compounds({cpd})
        r.set_left_compound_abundance(cpd, reactants[c])

    for i,c in enumerate(products):
        cpdname = modelseed.get_compound_by_name(c)
        if cpdname:
            cpd = CompoundWithLocation.from_compound(cpdname, 'c')
        else:
            PyFBA.log_and_message(f"biomass.py: No compound found for {c} in our compounds dataset", stderr=verbose)
            cpd = CompoundWithLocation(f"rctn{i}", c, 'c')
        r.add_right_compounds({cpd})
        r.set_right_compound_abundance(cpd, products[c])

    rcts = list(reactants.keys())
    prds = list(products.keys())
    rcts.sort()
    prds.sort()
    r.equation = " + ".join(["(" + str(reactants[x]) + ") " + x for x in rcts])
    r.equation += " > "
    r.equation += " + ".join(["(" + str(products[x]) + ") " + x for x in prds])

    return r
