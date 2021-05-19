# Parsers

We have parsers for a variety of different data formats and structures. Many of them feed directly into FBA pipelines
but some are accessory parsers that we have included. For example, we have parsers for the 
[Model SEED](http://www.theseed.org/models) data that should be installed in your [Biochemistry](../Biochemistry)
folder, and we have parsers for [SBML](http://www.sbml.org/) format files that you can bring in. In addition, in the
[media](../media) directory, we have defined a series of different media conditions and we have parsers to import those.

## Model SEED

The model SEED parser is broken into three main sections (although there are other components to the parser): 
*Compounds*, *Reactions*, and *Enzymes*. Each method returns a hash of the data stored with its `str()` function as the
key and the object as the value. This creates a convenient data structure for both looking up an object and iterating
over those objects.

The method `compounds_reactions_enzymes()` is the primary entry point to the model SEED data, and returns three data
structures, one each for `compounds`, `reactions`, and `enzymes`. These are the primary structures used to create and
analyse the FBA.

## SBML

The SBML parser uses the [beautiful soup](http://www.beatifulsoup.org/) XML parser to import all the data from an SBML
file. We have also provided code in [sbml_to_fba.py](../scripts/sbml_to_fba.py) to demonstrate how to extract the 
information from an SBML file, convert it to metabolism.Reaction and metabolism.Compound objects, and test for growth.

As an alternative, if you have `libsbml` installed you can also use the script 
[run_fba_sbml.py](../scripts/run_fba_sbml.py) to run the FBA from the SBML file using a different parser.


## Read Media

The [read_media.py](read_media.py) script parses a media file in a format that we designed and generates a list
of metabolism.Compound objects. The format is simple, tab separated text, with the columns [compound id, name, formula,
charge], with one compound per line.