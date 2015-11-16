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
important of which is the &Delta;G

If the reaction is truly reversible the probability can be 1 in both cases. If it is unidreictional the probability can
be 0 in one direction.

The likelihod that a reaction completes will be some product of its delta G and its p. We could also do something
simpler, e.g. if there  is a -ve delta G (favorable reaction) we can increase p and if there is a +ve delta G
(unfavorable reaction) we can decrease p.

We have arbitrarily called the reaction as proceeding from left to right but it could equally go the other way around.
The reaction.__eq__() method will check for either left->left/right->right and  left->right/right->left, and therefore
we include a reverse reaction method that will reverse a reaction.