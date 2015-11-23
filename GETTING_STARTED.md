# Getting started with PyFBA

The first thing you probably want to do is build a model for your genome. Because of the tight interplay between
RAST, SEED, and Model SEED, the easiest way to get started is to run your genome through RAST. 

## Using the SBML file

If you have done that and built the model using the Model SEED, you can download the SBML file from RAST, and 
try either [run_fba_sbml.py](example_code/run_fba_sbml.py) or [sbml_to_fba.py](example_code/sbml_to_fba.py). Both of 
these should give similar, but not identical answers to the answer that you got from the model 
SEED<sup id="a1">[1](#f1)</sup>.

## Using the genome annotation

If you have downloaded the annotation, there are two essential steps that you need to take to create a model:

1. Convert the genome annotation to reactions
2. Gapfill the reactions on different media.

### Convert the genome annotation to reactions

You can use the [example code](example_code) to get started. First, we will create a set of reactions from your genome:

```
fba_from_reactions.py 




---

<b id="f1">1</b> The reason that the answers are similar, but not identical, is because the linear solvers give
slightly different answers. The Model SEED uses a commercial linear solver, but you are probably using GLPK. [↩](#a1)