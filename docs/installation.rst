INSTALLING PyFBA
================

Prerequisites
=============

GLPK
----

To install PyFBA we first need to install a linear solver. We use the
`GNU Linear Programming Kit GLPK <https://www.gnu.org/software/glpk/>`__
program as our linear solver, although you can use others. In addition,
you will need to install the `Gnu Multiple
Precision <https://gmplib.org/>`__ arithmetic library. Finally, we
install PyGLPK, a Python wrapper around GLPK.

Install GLPK on Linux
~~~~~~~~~~~~~~~~~~~~~

Installing GLPK on Ubuntu or other Debian systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are on a Debian system you can install them with

``apt-get install libglpk-dev glpk-utils glpk-doc libglpk36 libgmp10 libgmpxx4ldbl``

(You maybe able to get away with less, but that is what I have
installed).

Once you have that installed, you should be able to continue with the
PyGLPK installation below

Installing GLPK on CentOS
^^^^^^^^^^^^^^^^^^^^^^^^^

If you are on a CentOS system, the easiest way to install this is to
first install GMP, download the GLPK library from GNU, and build it:

::

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

At this point, you need to add ``/usr/local/lib`` to your
``LD_LIBRARY_PATH`` and then continue with the PyGLPK installation below

``export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib python2.7 setup.py build python2.7 setup.py install``

Installing GLPK on Mac OSX 10.10 Yosemite
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unfortunately, the easiest way is to build glpk from source

Download the following file:
http://ftp.gnu.org/gnu/glpk/glpk-4.47.tar.gz and then open a Terminal
and install GLPK:

::

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

Installing GLPK on Windows
~~~~~~~~~~~~~~~~~~~~~~~~~~

You will need to install `WinGLPK <http://winglpk.sourceforge.net/>`__
but their installation instructions are not very clear. We will have
detailed installation instructions shortly.

Install PyGLPK
--------------

You need to install
`PyGLPK <https://github.com/bradfordboyle/pyglpk>`__. Make sure that you
install `this fork of
PyGLPK <https://github.com/bradfordboyle/pyglpk>`__ as the original has
not worked for quite some time!

Install the ModelSEEDDatabase
-----------------------------

We rely on the `Model SEED <http://www.theseed.org/models>`__ to provide
the biochemistry tables that we use (although we have designed PyFBA
with the potential to use biochemistry tables from other sources if they
become publicly available!). To get the latest ModelSEED database you
need to clone their GIT repository to a destination on your hard drive.
You then need to set the `ModelSEEDDatabase environment
variable <#set_the_environment_variables>`__ as explained below.

Python modules
--------------

PyFBA depends on a few different Python modules:

    * `libSBML <http://sbml.org/>`__
    * `Beautiful Soup 4 <http://www.crummy.com/software/BeautifulSoup/>`__
    * `PyGLPK <https://github.com/bradfordboyle/pyglpk>`__

As noted `above <#install_pyglpk>`__, you should install PyGLPK from
`GitHub <https://github.com/bradfordboyle/pyglpk>`__. However,
``setup.py`` will try and do the right thing for you.

libSBML and lxml
~~~~~~~~~~~~~~~~

One or two of the scripts (notably
`scripts/run\_fba\_sbml.py <scripts/run_fba_sbml.py>`__) require that
you have libSBML and lxml be installed to read the SBML files.

``setup.py`` will attempt to install these for you. If you wish to
install them manually you should be able to do so with ``pip install``:

::

        pip install python-libsbml-experimental
        pip install lxml

Beautiful Soup 4
~~~~~~~~~~~~~~~~

XML parsing is a pain in the butt, and so we use `Beautiful Soup
4 <http://www.crummy.com/software/BeautifulSoup>`__ to make life easy!
``setup.py`` should try and install this for you, but if you wish to do
it manually, you should be able to do so with ``pip intall``:

::

       pip install beautifulsoup4

Install PyFBA
=============

You should be able to install PyFBA from
`PyPI <https://pypi.python.org>`__ using ``pip install``:

::

        pip install pyfba

If that does not work, you can clone the git hub repository and run
setup.py manually:

::

        git clone https://github.com/linsalrob/PyFBA.git
        cd PyFBA
        # run the tests
        python setup.py test
        # install the code
        python setup.py install

If you do not have administrative (root) access to your machine you can
also install the code in a `local
directory <https://docs.python.org/2/install/#alternate-installation>`__:

::

        python setup.py install --user

Set the environment variables
=============================

PyFBA relies on the Model SEED Database, and we need to know where that
is installed. We also provide some example media files with the download
and you can set an environment variable (PYFBA\_MEDIA\_DIR) that points
to the location of those files if you want to include them for your
models.

Set the following two environment variables:

+-----------------------------+-----------------------------------------------------+-------------------------------------------------+
| Environment variable name   | points to                                           | example                                         |
+=============================+=====================================================+=================================================+
| ``ModelSEEDDatabase``       | the location of the Model SEED Database directory   | ``ModelSEEDDatabase=/data/ModelSEEDDatabase``   |
+-----------------------------+-----------------------------------------------------+-------------------------------------------------+
| ``PYFBA_MEDIA_DIR``         | the location of the media files                     | ``PYFBA_MEDIA_DIR=$HOME/FBA/media``             |
+-----------------------------+-----------------------------------------------------+-------------------------------------------------+

For more information on setting the environment variables, see one of
these sites: 

* `Windows <https://www.microsoft.com/resources/documentation/windows/xp/all/proddocs/en-us/sysdm_advancd_environmnt_addchange_variable.mspx>`__
* MacOS: 
    - You can set them globally as described `here <http://stackoverflow.com/questions/135688/setting-environment-variables-in-os-x/3756686#3756686>`__ or `here <https://developer.apple.com/library/mac/documentation/MacOSX/Conceptual/BPRuntimeConfig/Articles/EnvironmentVars.html>`__
    - You can set them in your .bashrc/.profile as described for Linux (below) 
    - You can set them in each terminal by using ``export ModelSEEDDatabase=/data/ModelSEEDDatabase`` 
* `Linux <http://www.cyberciti.biz/faq/set-environment-variable-linux/>`__

Tests
=====

The code in tests/testlp.py uses the example in the documentation (which
is also the example in the GLPK documentation) to solve a linear
programming problem, and only requires GLPX/GMP and PyGLPK. It does not
require any PyFBA code to solve a simple set of equations. You should
check that runs with ``nosetests tests/testlp.py`` and it should run a
single test that should pass. If that test does not pass, there is an
issue with your installation of GLPK, GMP, or PyGLPK and you should
check that each of them are installed in the correct locations.

There are many more tests in the `tests <PyFBA/tests/>`__ folder and you
can run all of them with ``nosetests tests/``. They should all run
without an error, and will test different aspects of the PyFBA
installation.

If you download and install the code from github, you can also run:

::

    python setup.py test

to run all the tests
