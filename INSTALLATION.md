# INSTALLING PyFBA

### Please Note

These installation instructions are for people who want to develop `PyFBA` and/or for some reason do not want to use
`conda`. 

If you want to use `PyFBA` for flux balance analysis, we strongly recommend you use `conda` to install it:

```commandline
conda create -n pyfba -c bioconda pyfba
conda activate pyfba
pyfba -v
pyfba help
```

# Installing for developers

If you want to install `PyFBA` for development, we recommend a simple process:
1. Use `conda` to install glpk. It really helps!
2. Clone this git repository 
3. Create a virtual environment.
4. Install the prerequisites
5. Build and install PyFBA.

```commandline
conda create -n glpk -c conda-forge glpk
conda activate glpk
git clone https://github.com/linsalrob/PyFBA.git
cd PyFBA
virtualenv venv
source venv/bin/activate
pip install -r requirements.txt
python setup.py install
```

This will create a virtural environment called `venv` and install the requirements and `PyFBA` in that environment.

_Note:_ I usually run `which pip` or `which python` after sourcing the virtual environment to make sure I am working
where I think I am.



# Installing without using conda

You don't _have_ to use conda, even though we highly, highly recommend it. The main problem will be installing `glpk`

## Prerequisites

### GLPK

To install PyFBA we first need to install a linear solver. We use the [GNU Linear Programming Kit
GLPK](https://www.gnu.org/software/glpk/) program as our linear solver, although you can use others. 

There are plenty of websites detailing how to install it. Older versions of this document detail 
installing it on CentOS, MacOS, and Windows. But now we use conda, and we have not maintained those instructions.

Once you have glpk installed, the instructions above should work.


# Tests

We have provided a test suite in [tests](PyFBA/tests/) folder and you can run all of them with `nosetests tests/`. They 
should all run without an error, and will test different aspects of the PyFBA installation.

If you download and install the code from GitHub, you can also run:

```
python setup.py test
```
to run all the tests
