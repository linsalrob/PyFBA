.. PyFBA documentation master file, created by
   sphinx-quickstart on Thu Nov 26 17:03:23 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyFBA
================

A python interface for Flux Balance Analysis

About PyFBA
+++++++++++

PyFBA is a Python flux-balance-analysis package that allows you to build models from genomes, gapfill models, and run
flux-balance-analysis on that model. The aim of PyFBA is to provide an extensible, Python-based platform for your
FBA work.

PyFBA is being developed by Daniel Cuevas, Taylor O'Connell, and Rob Edwards in `Rob's bioinformatics
group`_ at `San Diego State University`_ together with help from
Janaka Edirisinghe, Chris Henry, Ross Overbeek and others at `Argonne National Labs`_.

Installing PyFBA
++++++++++++++++

To use PyFBA you need Python 2.7 or greater, and you need to install the GNU GLPK and a Python
wrapper for that program, `pyGLPK`_ available from github. See the `installation`_ page for more details.

We also leverage the `Model SEED`_ repository to get all the latest biochemistry tables. You should install that
somewhere on your machine, and set a `ModelSEEDDatabase environment` variable that points to that directory so we know
where you have installed it.

Our `installation`_ page has detailed instructions on installing PyFBA and getting everything running.

.. _Rob's bioinformatics group: http://edwards.sdsu.edu/research
.. _San Diego State University: http://www.sdsu.edu
.. _Argonne National Labs: http://www.theseed.org/
.. _installation: installation.html
.. _pyGLPK: https://github.com/bradfordboyle/pyglpk
.. _Model SEED: https://github.com/ModelSEED/ModelSEEDDatabase.git


Getting Started with PyFBA
++++++++++++++++++++++++++

Once you have installed GLPK, PyGLPK, and PyFBA, you will most likely want to build a model from a genome, gap fill that
model, and test it for growth on different media. We have detailed instructions on `getting started`_ to that walk you
through the step-by-step procedures that you need to use to run flux balance analysis on your own genome.

.. _getting started: getting_started.html


Explore the API
+++++++++++++++

We have developed an extensive API that is centered around two concepts:
    * Getting your FBA up and running easily
    * Allowing you to extend, update, and adapt the API to meet your requirements

To explore the API, check out the details:

.. toctree::
   :maxdepth: 2
   
   api.rst



You can also find what you are looking for directly:
++++++++++++++++++++++++++++++++++++++++++++++++++++

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

