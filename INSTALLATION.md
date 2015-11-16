# INSTALLING PyFBA

## GLPK

To install PyFBA we first need to install a linear solver. We use the [GNU Linear Programming Kit
GLPK](https://www.gnu.org/software/glpk/) program as our linear solver, although you can use others. In addition, you
will need to install the [Gnu Multiple Precision](https://gmplib.org/) arithmetic library. Finally, we install PyGLPK, a
Python wrapper around GLPK.

### Install GLPK on Linux

#### Installing GLPK on Ubuntu or other Debian systems

If you are on a Debian system you can install them with

``` apt-get install libglpk-dev glpk-utils glpk-doc libglpk36 libgmp10 libgmpxx4ldbl ```

(You maybe able to get away with less, but that is what I have installed).

Once you have that installed, you should be able to continue with the PyGLPK installation below

#### Installing GLPK on CentOS

If you are on a CentOS system, the easiest way to install this is to first install GMP, download the GLPK library from
GNU, and build it:

    yum install gmp.x86_64 gmp-devel.x86_64
    mkdir glpk 
    cd glpk/ 
    wget ftp://mirrors.kernel.org/gnu/glpk/glpk-4.56.tar.gz
    tar zxf glpk-4.56.tar.gz
    cd glpk-4.56 
    ./configure  --with-gmp
    make 
    make test 
    make install

At this point, you need to add `/usr/local/lib` to your `LD_LIBRARY_PATH` and then continue with the PyGLPK installation
below

``` export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib python2.7 setup.py build python2.7 setup.py install ```


### Installing GLPK on Mac OSX 10.10 Yosemite

Unfortunately, the easiest way is to build glpk from source

Download the following file: http://ftp.gnu.org/gnu/glpk/glpk-4.47.tar.gz and then open a Terminal and install GLPK:

    mkdir glpk_install
    cd glpk_install
    # Move tar file here to this directory 
    mv ~/Downloads/glpk-4.47.tar.gz ./ 
    tar xzf glpk-4.47.tar.gz 
    cd glpk-4.47
    ./configure
    make 
    make test 
    make install

### Installing GLPK on Windows

You will need to install [WinGLPK](http://winglpk.sourceforge.net/) but their installation instructions are not very
clear. We will have detailed installation instructions shortly.

## Install PyGLPK

You need to install [PyGLPK](https://github.com/bradfordboyle/pyglpk). Make sure that you install [this fork of
PyGLPK](https://github.com/bradfordboyle/pyglpk) as the original has not worked for quite some time!


## Install PyFBA

To install PyFBA use the `git clone` link you see on the right hand side:
 
```
    git clone https://github.com/linsalrob/PyFBA.git
```
 
This will download the PyFBA code and install it in a directory called (imaginatively) **PyFBA**. 

We rely on the [Model SEED](http://www.theseed.org/models) to provide the biochemistry tables that we use (although we 
have designed PyFBA with the potential to use biochemistry tables from other sources if they become publicly 
available!). To get the latest ModelSEED data, you need to run these two `git` commands in your new `PyFBA` directory:

```
    git submodule init
    git submodule update
```

These will checkout and download all the Model SEED data for you.

As an alternative to this two step process, you can also just use the one line recursive git clone:

```
    git clone --recursive https://github.com/linsalrob/PyFBA.git
```

which will download all of the PyFBA *and* Model SEED code in one shot.



# Tests

The code in tests/testlp.py uses the example in the documentation (which is also the example in the GLPK documentation)
to solve a linear programming problem, and only requires GLPX/GMP and PyGLPK. It does not require any PyFBA code to
solve a simple set of equations. You should check that runs with `nosetests tests/testlp.py` and it should run a single
test that should pass. If that test does not pass, there is an issue with your installation of GLPK, GMP, or PyGLPK and
you should check that each of them are installed in the correct locations.

There are many more tests in the [tests](tests/) folder and you can run all of them with `nosetests tests/`. They 
should all run without an error, and will test different aspects of the PyFBA installation.

