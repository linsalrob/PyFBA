import sys
import PyFBA
from PyFBA import lp, log_and_message


def create_stoichiometric_matrix(reactions_to_run, reactions, compounds, media, biomass_equation,
                                 uptake_secretion=None, verbose=False):
    """Given the reactions data and a list of RIDs to include, build a
    stoichiometric matrix and load that into the linear solver.

    The uptake and secretion reactions (sometimes called boundary reactions) are reactions that allow
    compounds to flow into the media and secreted compounds away from the cell. You can either provide these (e.g. if
    you are parsing an XML (SBML) file, or we will calculate them for you based on your media and reactions.

    We also take this opportunity to set the objective function (as it is a member of the SM).

    :param uptake_secretion: An optional hash of uptake and secretion reactions that should be added to the model
    :type uptake_secretion: a dict of Reactions
    :param compounds: a set of the compounds present
    :type compounds: set[PyFBA.metabolism.CompoundWithLocation]
    :param reactions: a dict of all the reactions
    :type reactions: dict[PyFBA.metabolism.Reaction]
    :param reactions_to_run: just the reaction ids that we want to include in our model
    :type reactions_to_run: set
    :param media: a set of compounds that would make up the media
    :type media: set
    :param biomass_equation: the biomass_equation equation as a Reaction object
    :type biomass_equation: metabolism.Reaction
    :param verbose: print more information
    :type verbose: bool
    :returns: Sorted lists of all the compounds and reactions in the model, and a revised reactions dict that includes the uptake and secretion reactions
    :rtype: list, list, dict

    """

    sm = {}  # the matrix
    reaction_cpds = set()  # all the cpds in the matrix

    # unify compounds from dicts and sets
    if type(compounds) == dict:
        log_and_message("create_stoichiometric_matrix is converting compounds from a dict to a set")
        compounds = set(compounds.values())

    # initialize the compounds set and the stoichiometric matrix with
    # everything in the media. The order of compounds is irrelevant
    for c in media:
        if c not in compounds:
            log_and_message(f"create_stoichiometric_matrix is adding compound {c} from media to compounds", stderr=True)
            compounds.add(c)
        if verbose and type(c) != PyFBA.metabolism.compound.CompoundWithLocation:
            log_and_message(f"create_stoichiometric_matrix is parsing the media, {c} is a {type(c)} " +
                            f"(not a cpod with location)", stderr=True)
        reaction_cpds.add(c)
        sm[c] = {}

    # iterate through the reactions
    for r in reactions_to_run:
        for c in reactions[r].left_compounds:
            reaction_cpds.add(c)
            if c not in sm:
                sm[c] = {}
            if not isinstance(c, PyFBA.metabolism.compound.CompoundWithLocation):
                log_and_message(f"Warning. In parsing left compounds for the SM, {c} is a {type(c)}", stderr=verbose)
            sm[c][r] = 0 - reactions[r].get_left_compound_abundance(c)

        for c in reactions[r].right_compounds:
            reaction_cpds.add(c)
            if c not in sm:
                sm[c] = {}
            if not isinstance(c, PyFBA.metabolism.compound.CompoundWithLocation):
                log_and_message(f"Warning. In parsing right compounds for the SM, {c} is a {type(c)}", stderr=verbose)
            sm[c][r] = reactions[r].get_right_compound_abundance(c)

    for c in biomass_equation.left_compounds:
        if c not in compounds:
            compounds.add(c)
            log_and_message(f"create_stoichiometric_matrix added compound {c} from biomass to compounds", stderr=True)
        if verbose and type(c) != PyFBA.metabolism.compound.CompoundWithLocation:
            log_and_message(f"In parsing biomass left, {c} is a {type(c)}", stderr=True)
        reaction_cpds.add(c)
        if c not in sm:
            sm[c] = {}
        sm[c]["BIOMASS_EQN"] = 0 - biomass_equation.get_left_compound_abundance(c)
    for c in biomass_equation.right_compounds:
        if c not in compounds:
            compounds.add(c)
            log_and_message(f"create_stoichiometric_matrix added compound {c} from biomass to compounds", stderr=True)
        if verbose and type(c) != PyFBA.metabolism.compound.CompoundWithLocation:
            log_and_message(f"In parsing biomass right, {c} is a {type(c)}", stderr=True)
        reaction_cpds.add(c)
        if c not in sm:
            sm[c] = {}
        sm[c]["BIOMASS_EQN"] = biomass_equation.get_right_compound_abundance(c)

    # Add the uptake/secretion reactions. These are reactions that allow things to flow from the media
    # into the reaction, or from the cell outwards.
    #
    # The reactions are determined by the external compounds (and biomass_equation), but we only add the left side
    # of the equation to the stoichiometric matrix. This means that they can appear and disappear at will!
    #
    # When we set the reaction bounds we determine which things are in the media unless they are provided for you

    if not uptake_secretion:
        uptake_secretion = PyFBA.fba.uptake_and_secretion_reactions(reaction_cpds)
        log_and_message(f"create_stoichiometric_matrix found {len(uptake_secretion)} uptake and secretion reactions",
                        stderr=True)
    for r in uptake_secretion:
        reactions[uptake_secretion[r].id] = uptake_secretion[r]
        for c in uptake_secretion[r].left_compounds:
            reaction_cpds.add(c)
            if c not in sm:
                sm[c] = {}
            sm[c][uptake_secretion[r].id] = 0 - uptake_secretion[r].get_left_compound_abundance(c)

    # now we need to make this into a matrix sorted by
    # reaction id and by cpds
    cp = list(reaction_cpds)
    cp.sort()
    rc = list(reactions_to_run)
    rc.sort()
    rc += [uptake_secretion[x].id for x in uptake_secretion]

    # it is important that we add these at the end
    rc.append("BIOMASS_EQN")

    if verbose:
        sys.stderr.write(f"In the model there are : {len(cp)} compounds and {len(rc)} reactions\n")

    # here we create the matrix from our sm hash
    data = []
    for i, j in enumerate(cp):
        if j not in sm:
            sys.exit("Error while parsing: no " + j + " in sm")
        data.append([])
        for c in rc:
            if c in sm[j]:
                data[i].append(sm[j][c])
            else:
                data[i].append(0.0)

    # load the data into the model
    PyFBA.lp.load(data, [str(c) for c in cp], [str(r) for r in rc], verbose=verbose)

    # now set the objective function.It is the biomass_equation
    # equation which is the last reaction in the network
    ob = [0.0 for r in rc]
    ob[-1] = 1

    PyFBA.lp.objective_coefficients(ob)

    return cp, rc, reactions
