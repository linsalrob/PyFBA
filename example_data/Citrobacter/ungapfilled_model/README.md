# Raw Data

The data in this directory is raw data from *Citrobacter sedlakii* that you can use to test some of the 
[gap-filling](../../../gapfill) approaches that we describe:

* [citrobacter.assigned_functions](citrobacter.assigned_functions) is a list of functional roles from RAST ID 67826.8
* [citrobacter.reactions](citrobacter.reactions] is a list of reactions the assigned functions using roles_to_reactions
* [citrobacter.roles](citrobacter.roles) is a list of roles in *all other Citrobacter* that can be used for gap-filling
* [closest.genomes](closest.genomes) is a list of closest genomes from RAST
* [closest.genomes.roles](closest.genomes.roles) is a list of all the roles in the closest genomes that can be used 
for gapfilling.

To gap fill this model we will use a wide range of techniques. `pyfba gapfill_roles` uses the list of reactions in 
citrobacter.reactions to start the initial model. We then propose reactions from a variety of different sources, 
including both the closest.genomes.roles file and the citrobacter.roles file. After successfully getting growth, we
start pruning the reactions using the methods in the [gapfilling bisection](../../../gapfill/bisection.py) library
that trims out unwanted reactions.

To run this gap filling approach, us the command:

```
    python gapfill_roles -r example_data/Citrobacter/ungapfilled_model/citrobacter.reactions \
    -m MOPS_NoC_Alpha-D-Glucose.txt -c example_data/Citrobacter/ungapfilled_model/closest.genomes.roles \
    -g example_data/Citrobacter/ungapfilled_model/citrobacter.roles
```
