# Metaboolism

The metabolism objects are central to our model and consist of three main classes:

* [Reaction](reaction.py)
* [Compound](compound.py) 
* [Enzyme](enzyme.py)

## Reaction

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


    The left and right compounds are sets, which should contain compound
    objects. 

    The left and right abundances are hashes, with the key being the
    compound and the value being the amount of that compound produced or
    consumed by the reaction.
    
    The two probabilities are from left to right (pLR) and right to left
    (pRL). These are relative to how the compounds are defined.

    The set of enzyme(s) that complete the reaction should be included.
    Usually (sometimes? often?) this will be a single enzyme but 
    sometimes it may be several enzymes. (I may refactor this to be an
    enzyme object).

    Pegs are a set of proteins that are involved in fullfulling this 
    reaction.

    Input (inp) and output (outp) are the reactions that will start 
    everything off. Generally, import reactions (converting 
    extracellular compounds to intracellular compounds) are inputs. 
    Converting intracellular compounds to extracellular are output

    ran is a boolean flag you can set to note if this reaction
    was ran, so you can then iteration through all reactions
    to see which ones we ran. Feel free to set/unset this
    flag at will.

    is_biomass_reaction and biomass_direction indicate whether this
    is a biomass reaction (T/F) and if so the direction that it runs
    (L->R or R->L)

    Is_gapfilled indicates whether this reaction was added by gapfilling










We have arbitrarily called the reaction as proceeding from left to right but it could equally go the other way around.
The reaction.__eq__() method will check for either left->left/right->right and  left->right/right->left, and therefore
we include a reverse reaction method that will reverse a reaction.






<b id="f1">1</b> I realize that these are not Pythonic names, but plr and prl just seem wrong when talking about the
probability of moving left to right and right to left. [↩](#a1)