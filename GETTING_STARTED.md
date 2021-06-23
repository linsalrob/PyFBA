# Getting started with PyFBA

The first thing you probably want to do is build a model for your genome. Because of the tight interplay between
PyFBA, PATRIC, RAST, SEED, and Model SEED, the easiest way to get started is to run your genome through RAST or PATRIC.  

We are going to assume that you have annotated your genome through the [PATRIC annotation 
pipeline](https://www.patricbrc.org/). Once you have your annotated genome, download the `features` file that will
be called `_<GENOME NAME>_.features.txt`

This file should have five columns:
feature_id	location	type	function	aliases	protein_md5

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

You can use the [example code](example_code) to get started. First, we will create a set of reactions from your genome.
First, download the assigned_functions file from the genome directory, or get a list of all funcational roles in your
genome. Next, convert those to a list of reactions, using one of two commands:

If you have an assigned_functions file that has [peg, functional role]:

```
    python example_code/assigned_functions_to_reactions.py -a assigned_functions > reactions.list
```

or, if you just have a list of functional roles, one per line:

```
    python example_code/assigned_functions_to_reactions.py -r functional_roles > reactions.list
```

### Build a gapfilled model

We build a model from the reactions and then try and gap fill it using all of our approaches. 

```
python scripts/gapfill_from_reactions.py -r reactions.list -m MOPS_NoC_Alpha-D-Glucose.txt > out.txt 2> out.err
```

This will use the media file `MOPS_NoC_Alpha-D-Glucose.txt` and try and gap fill the model so that it grows on this 
media.

---


<b id="f1">1</b> The reason that the answers are similar, but not identical, is because the linear solvers give
slightly different answers. The Model SEED uses a commercial linear solver, but you are probably using GLPK. [↩](#a1)