# Gap-filling

## Overview

Gap-filling is an essential process in FBA. The basic model that we build typically does not have enough reactions to
grow. Therefore, we need to identify the reactions that are missing from the model, add them to the model, and test 
for growth.

There are a number of different ways that you can propose reactions to be added to the model, and most of the code in 
this package is associated with making suggestions for different reactions that could be added to a model and that 
could result in growth.

The steps that we take are
    
    1. Identify a list of reactions that could be added to the model
    2. Check which reactions are already in the model, and exclude them
    3. Propose a set of new reactions that could be added
    4. Test whether the model with the new reactions results in growth.
    
The methods here, therefore, almost all return a set of reactions that are not currently in the model but which should
be added to the model.

In the standard [Model SEED Biochemistry files](https://github.com/ModelSEED/ModelSEEDDatabase) there over 30,000 reactions. This
is too many reactions that we can just add them all and then test which ones should be there (i.e. a complete brute
force approach), so we finesse it by adding just those reactions that we have evidence for being associated with the
model.

We can break those 34,000 reactions into two groups:

    * The set of reactions for which we have proteins known to perform the reaction
    * The set of reactions for which there is no known protein (or we have not annotated it)

The first set is approximately 2,000 reactions, and adding these often involves using homology approaches to suggest
the reactions that should be added to your genome.

The second set of reactions is about 32,000 reactions, and arise for different reasons. Some reactions are spontaneous
and do not require proteins to catalyze the reactions. Some reactions are known and proven, but we don't have the
protein identified that performs the reaction. Some reactions are hypothesized to happen because biochemistry and
thermodynamics says that they should, but they have not been shown yet. Regardless for why reactions are in this set
we need mechanisms to filter the ~32,000 reactions to a manageable set that we can test.

# Gap-filling methods

The currently implemented methods in gap-filling are:

### Essential reactions

There are a set of [essential reactions](essentials.py) that don't have proteins associated with them. From our 
analysis of about 2,000 gapfilled models, all of the models have these reactions! Therefore, we just include these 
as suggestions to ensure that your model has them.

### Media import reactions

If your model doesn't grow on the media, and you know it should, one potential problem could be that we have not
identified the transporters. Transporter proteins are notoriously easy to identify but almost impossible to determine
the correct compound that is being transported just based on sequence homology. (Sometimes, however, you can identify
the substrate from genes that are adjacent to the transporter.) The [media reactions](media.py) suggests reactions that
ensure import of substrates in your media formulation that are not already being transported by other routes.

### Orphan compounds

Compounds can not be consumed if they are not produced, and most compounds that are produced are consumed (or excreted)
we therefore [check for compounds](orphan_compounds.py) that are only present in single reactions (orphan compounds),
and suggest reactions to add to your model based on these compounds.

### Probability of inclusion

The probability is based off the idea that the larger a percentage of compounds that a reaction has present in your 
model, the more likely that reaction is to be present in your model. We [propose reactions to add](probability.py)
based on the compounds that you already have and how the proposed reactions would fit in.

### Other genomes and organisms

Because of the way evolution works, your organism is likely to be similar to other organisms nearby on the phylogenetic 
tree. Therefore, a successful approach to gap-filling also involves identifying the roles that are present in 
similar<sup id="a1">[1](#f1)</sup> organisms. We start with a file that has tab-separated file that has two columns:

    * role
    * fraction of organisms that have that role

and based on the roles in the first column (optionally filtered by their frequency of occurrence in the second column),
we [propose reactions to add to your model](roles.py).

### Subsystems

The reactions in your model map to subsystems, but not all of those subsystems will be complete. Therefore, we [propose
reactions to add](subsystems.py) to your model that will complete the subsystems. You can optionally provide a 
threshold to limit to only those reactions that have a certain level of completeness.

### Reactions that [do or do not] map to proteins

As a last resort, we can add all reactions that either map to proteins, or all reactions that do not map to proteins. 
This should result in growth, but you may end up with such a large model that the computation because too expensive.

If you try either of these approaches it is critical that you prune the reactions as we have described below.

# Pruning the suggested reactions

There are many different ways to suggest reactions that should be added, and a lot of potential reactions that could
be added, but it clearly doesn't make sense to add reactions for which none of the compounds are currently in the model.
This would just add spurious reactions. We therefore have [code to limit the reactions](limit_reactions.py) to just
those reactions that include compounds already in your model.

Once we have some suggestions, we need to refine the set of suggestions to just those reactions that are
required for growth. the [bisections](bisections.py) code has methods for taking lists of reactions and splitting them
into different parts (either equally or unequally) and returning two separate sets of reactions. We also try and 
optimize the splits based on clustering reactions, for example by the compounds that are present in the reaction.


# Writing new gapfilling methods

Writing new methods is trivial. The only constraint is that you should return a set of reaction IDs that are not 
currently in the model. 

A typical gap filling approach will be to accept as parameters to your methods the reactions dict that has all reaction
data, and the set of reactions_to_run that are currently in the model (so you know which reactions are unique). Identify
which reactions you think should be added to the model, and return that set of IDs.


---

<b id="f1">1</b> Your definition of similar may vary, but we often use all organisms in the same genus or all
closely related organisms, e.g. as suggested by [RAST](http://rast.nmpdr.org/). [↩](#a1)