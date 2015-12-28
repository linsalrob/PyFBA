

def suggest_essential_reactions():
    """
    Identify a set of reactions that you should add to your model for growth because they are essential reactions.

    There are 110 reactions (one of which is biomass) that are in _every_ model produced thus far and we include them
    in all models. This is a set of those reactions just to make sure that we have added them!

    :return: A set of proposed reactions that should be added to your model to see if it grows
    :rtype: set
    """

    rxns = {
        "rxn13784", "rxn13783", "rxn13782", "rxn12224", "rxn12008", "rxn11946", "rxn10571", "rxn10473", "rxn10338",
        "rxn10337", "rxn10336", "rxn10266", "rxn10265", "rxn10260", "rxn10259", "rxn10233", "rxn10232", "rxn10227",
        "rxn10226", "rxn10221", "rxn10220", "rxn10215", "rxn10214", "rxn10206", "rxn10205", "rxn10199", "rxn08333",
        "rxn05667", "rxn05651", "rxn05555", "rxn05468", "rxn05467", "rxn05454", "rxn05452", "rxn05406", "rxn05405",
        "rxn05404", "rxn05402", "rxn05401", "rxn05400", "rxn05398", "rxn05397", "rxn05396", "rxn05394", "rxn05393",
        "rxn05392", "rxn05390", "rxn05389", "rxn05388", "rxn05386", "rxn05385", "rxn05384", "rxn05383", "rxn05381",
        "rxn05380", "rxn05379", "rxn05377", "rxn05376", "rxn05375", "rxn05373", "rxn05372", "rxn05371", "rxn05369",
        "rxn05368", "rxn05367", "rxn05365", "rxn05364", "rxn05363", "rxn05361", "rxn05360", "rxn05359", "rxn05358",
        "rxn05319", "rxn05195", "rxn05116", "rxn05064", "rxn05029", "rxn04457", "rxn04456", "rxn04139", "rxn04133",
        "rxn04132", "rxn03904", "rxn03901", "rxn03893", "rxn03538", "rxn03537", "rxn03408", "rxn03397", "rxn03395",
        "rxn03393", "rxn03164", "rxn03150", "rxn03012", "rxn02916", "rxn02897", "rxn02666", "rxn02374", "rxn02286",
        "rxn02285", "rxn02056", "rxn02011", "rxn02008", "rxn01664", "rxn01208", "rxn00851", "rxn00461", "rxn00392",
        "rxn00062"
    }

    return rxns
