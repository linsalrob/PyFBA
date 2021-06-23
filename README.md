[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.sdsu.edu/research)
[![DOI](https://www.zenodo.org/badge/46241465.svg)](https://www.zenodo.org/badge/latestdoi/46241465)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/PyFBA)
[![PyPi](https://img.shields.io/pypi/pyversions/pyfba.svg?style=flat-square&label=PyPi%20Versions)](https://pypi.org/project/PyFBA/)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/pyfba.svg?style=flat-square&label=BioConda%20install)](https://anaconda.org/bioconda/pyfba)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/pyfba/badges/version.svg)](https://anaconda.org/bioconda/pyfba)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/pyfba/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)


# PyFBA
A python implementation of flux balance analysis to model microbial metabolism. [Read more ...](http://linsalrob.github.io/PyFBA/)

## About PyFBA

PyFBA is a Python flux-balance-analysis package that allows you to build models from genomes, gapfill models, and run
flux-balance-analysis on that model. The aim of PyFBA is to provide an extensible, python-based platform for
FBA work.

PyFBA is being developed by Daniel Cuevas, Taylor O'Connell, and Rob Edwards in [Rob's bioinformatics
group](http://edwards.flinders.edu.au/) initially at [San Diego State University](http://www.sdsu.edu/) and now at 
[Flinders University](http://www.flinders.edu.au/). Amazing help was also provided by the developers of the ModelSEED, 
in particular  Janaka Edirisinghe, Chris Henry, Ross Overbeek and others at [Argonne National Labs](http://www.theseed.org/).

You can [read more about PyFBA](http://linsalrob.github.io/PyFBA/) on our github.io pages.

## Installing PyFBA

PyFBA is available in [bioconda](https://bioconda.github.io/) and we recommend installing it that way. It's what we 
use!

```commandline
conda create -n pyfba -c bioconda pyfba
conda activate pyfba
pyfba -v
pyfba help
```

There are other options described in the [installation documents](INSTALLATION.md), but just use conda!

## Getting Started with PyFBA

Once you have installed PyFBA, you will most likely want to build a model from a genome, gap fill that
model, and test it for growth on different media. We have [detailed instructions](GETTING_STARTED.md) that walk you through the step-by-step
procedures that you need to use to run flux balance analysis on your own genome.

### Citing PyFBA

Please use the command:

```commandline
pyfba citations
```

to get the citations for PyFBA. They are available in _plain text_ and _bibtex_ format, please let us know if you would like
other formats!

### Copyright and License

PyFBA is copyright Daniel Cuevas, Taylor O'Connell, and Rob Edwards, and is released under the MIT license.
