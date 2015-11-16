# Metaboolism

The metabolism objects are central to our model and consist of three main classes:

* [Reaction](#reaction)
* [Compound](#compound) 
* [Enzyme](#enzyme)

In addition, we have some accessory classes that live here and define some metabolic processes that we use:

* [biomass](#biomass)

## [Reaction](reaction.py)

The reaction class contains all the information about a specific reaction. A reaction is the central concept of metabolism and is
the conversion of substrates to products. All of the substrates and the products are [Compound](Compound.py) objects.

A reaction is an object that describes how to get from one compound to another. We need to know what the compound(s) on
the left of the equation are, what the compounds on the right of the reaction are, and the probability that the reaction
proceeds in either direction. The probability can come from a number of different measures or options, the most
important of which is the &Delta;G.

If the reaction is truly reversible the probability can be 1 in both cases. If it is unidreictional the probability can
be 0 in one direction, and we can easily test the probablity using the pLR and pRL variables<sup id="a1">[1](#f1)</sup>.
In the future we intend that the likelihood that a reaction completes will be some product of its &Delta;G and its p,
however at the moment we have not implemented this capability.

At a bare minimum we need a name for the reaction. The name can either be the reaction id (e.g. modelSEED or KEGG id), 
or another name for this reaction. We tend to use the ID of the [biochemistry](../Biochemistry) that we are parsing,
and thus currently use the Model SEED ID.

The equation is the full formula of the reaction. It is optional but you should include it if you have it because it is
used in some of the output options.

The direction is the direction that the equation can run. Acceptable values are:

```
        None    We don't know the direction
        >       Left to right
        <       Right to left
        =       Bidirectional
```


The left and right compounds are sets of [Compound](compound.py) objects.

The left and right abundances are hashes, with the key being the compound and the value being the amount of that
compound produced or consumed by the reaction.

The two probabilities are from left to right (pLR) and right to left (pRL). These are relative to how the compounds are
defined.

The complexes that make up the enzymes that perform the reaction should be included so that we know which reaction goes
with which complex.

Pegs are a set of proteins that are involved in fullfulling this reaction.

Input (inp) and output (outp) are the reactions that will start everything off. Generally, import reactions (converting
extracellular compounds to intracellular compounds) are inputs. Generally, converting intracellular compounds to
extracellular are output. This can be computed using `check_input_output()` which will test whether compounds are
imported or exported.

is_biomass_reaction and biomass_direction indicate whether this is a biomass reaction (T/F) and if so the direction that
it runs (L->R or R->L)

Is_gapfilled indicates whether this reaction was added by gapfilling, and gapfill_method allows us to note which
gapfilling method was used to identify the reaction.

The is_uptake_secretion flag is used to determine whether this is an uptake/secretion reaction that is provided to
ensure that we have enough media components.

We have arbitrarily called the reaction as proceeding from left to right but it could equally go the other way around.
The reaction.__eq__() method will check for either left->left/right->right and  left->right/right->left, and therefore
we include a reverse reaction method that will reverse a reaction.


## [Compound](compound.py)

A compound is a metabolite in our model, and is represented by a name and a location. Note that we typically use the
name and not the compound ID so that we can easily look at the compound! We provide a model_seed_id variable that
can be set with alternate IDs as required, and an abbreviation that can be used to store different names.

The location is generally one of

```
    e: extracellular
    c: cytoplasmic
    h: chloroplast
    p: periplasm
```

for most of our reactions they tend to focus on extracellular or cytoplasmic.


Compounds are involved in reactions, and we have a set of reactions that these compounds are connected to.

The formula is the normal chemical formula for the compound, and the charge associated wtih the compound can also be 
provided. The `calculate_molcular_weight()` method has not yet been implemented but is on the todo list. 

We have two booleans, `common` denotes common compounds that are present in many reactions. Often we want to parse
through lists of compounds and the reactions that they are involved in, but we don't want to parse through, for example
H<sub>2</sub>O which is present in almost every reaction! The `common` flag denotes compounds that are in 
plenty<sup id="a2">[2](#f2)</sup>. of reactions and should be skipped.

The other boolean is `uptake_secretion` that denotes whether the compound is involved in uptake from the media or
secretion back to the media.


## [Enzyme](enzyme.py)

The Enzymes class connects [Reaction](reaction.py) objects with functional roles. In the Model SEED  a set of 
functional roles contribute to a complex in a many to many relationship (many functional roles can be in one complex and
one functional role can be in many complexes). Similarly, complexes are connected to reactions in a many to many 
relationship: many complexes can perform one reaction and one complex can be in many reactions. The complex is thus the
hub that joins reactions and functional roles.

When we build a model from a genome we start with functional roles. This means that we have to convert those to 
complexes and from there extract the reactions that should be in our model. We do this by way of the enzyme objects.

The enzyme is not a complicated object, just a connector between reactions and roles, and you can think about is as 
comprising: The subunit(s) that make up the enzyme, the genes that encode those subunit(s), and  the reactions that
this enzyme is connected to.


### [Biomass](biomass.py)

The biomass reactions are found in the [Model SEED](Biochmistry/ModelSEEDDatabase/Templates) database, but we have 
abstracted them here so that we can manipulate them (not all compounds are required in all biomass equations), and 
to demonstrate how you should go about creating a biomass object for the FBA.

---



<b id="f1">1</b> I realize that these are not Pythonic names, but plr and prl just seem wrong when talking about the
probability of moving left to right and right to left. [↩](#a1)

<b id="f2">2</b> Plenty is a user-definable number which is currently set at 50. [↩](#a2)