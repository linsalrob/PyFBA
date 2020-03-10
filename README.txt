

PYFBA


A python implementation of flux balance analysis to model microbial
metabolism. Read more ...


About PyFBA

PyFBA is a Python flux-balance-analysis package that allows you to build
models from genomes, gapfill models, and run
flux-balance-analysis on that model. The aim of PyFBA is to provide an
extensible, python-based platform for
FBA work.

PyFBA is being developed by Daniel Cuevas, Taylor O'Connell, Shane Levi,
and Rob Edwards in Rob's bioinformatics
group at San Diego State University together with help from
Janaka Edirisinghe, Chris Henry, Ross Overbeek and others at Argonne
National Labs.

You can read more about PyFBA on our github.io pages.


Installing PyFBA

To use PyFBA you need Python 3.0 or greater, and you need to install the
GNU GLPK and a Python
wrapper for that program, pyGLPK available from github.

We also leverage the Model SEED GitHub repository with all the
latest biochemistry tables.

Our installation page has detailed instructions on installing PyFBA and
getting everything running.


Getting Started with PyFBA

Once you have installed GLPK, PyGLPK, and PyFBA, you will most likely
want to build a model from a genome, gap fill that
model, and test it for growth on different media. We have detailed
instructions that walk you through the step-by-step
procedures that you need to use to run flux balance analysis on your own
genome.


Updating to Python 3

This branch of PyFBA is compatible with Python 3, however that currently
breaks compatibility with Python 2 because we use
an open syntax that is not supported in earlier Python versions.

To update to Python 3 from a working Python 2 installation, you will
need to reinstall pyGLPK available from
github. Checkout that code, cd into the directory, and then install it
using python3 setup.py. Everything else should work just fine.

If you cd into the PyBFA directory, nosetests3 tests/ should
successfully run 77 tests with no errors.

Copyright and License

PyFBA is copyright Daniel Cuevas, Taylor O'Connell, and Rob Edwards, and
is released under the MIT license.
