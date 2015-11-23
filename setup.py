import os
from distutils.core import setup
from setuptools.command.test import test as TestCommand

import sys

import PyFBA
import io

# this is taken from https://www.jeffknupp.com/blog/2013/08/16/open-sourcing-a-python-project-the-right-way/ and
# allows README.txt and CHANGES.txt to automatically be used by PyPI. Create README.txt using http://pandoc.org/:
# pandoc -t plain -o README.txt README.md

if 'ModelSEEDDatabase' not in os.environ:
    sys.stderr.write("Please ensure that you install the Model SEED Database somewhere, and set the environment " +
                     "variable ModelSEEDDatabase to point to that directory.\n" +
                     " See INSTALLATION.md for more information\n")
    sys.exit(-1)

if not os.path.exists(os.environ['ModelSEEDDatabase']):
    sys.stderr.write("The ModelSEEDDatabase environment variable points to {}".format(os.environ['ModelSEEDDatabase']) +
                     " but that location does not exist\n")
    sys.exit(-1)

if not os.path.exists(os.path.join(os.environ['ModelSEEDDatabase'], 'Biochemistry/reactions.master.tsv')):
    sys.stderr.write("The ModelSEEDDatabase at {} ".format(os.environ['ModelSEEDDatabase']) +
                     "does not seem to be complete. Did you check it out from GitHub?\n")
    sys.exit(-1)


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.txt', 'CHANGES.txt')

# Allow nosetests for our python setup.py test framework

# Inspired by the example at https://pytest.org/latest/goodpractises.html
class NoseTestCommand(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        # Run nose ensuring that argv simulates running nosetests directly
        import nose
        nose.run_exit(argv=['tests'])

setup(
    name='PyFBA',
    version=PyFBA.__version__,
    packages=['PyFBA', 'PyFBA.lp', 'PyFBA.fba', 'PyFBA.parse', 'PyFBA.tests', 'PyFBA.filters', 'PyFBA.gapfill',
              'PyFBA.metabolism'],
    url='http://linsalrob.github.io/PyFBA/',
    license='The MIT License (MIT)',
    author='Rob Edwards',
    author_email='raedwards@gmail.com',
    long_description=long_description,
    platforms='any',
    test_suite = 'nose.collector',
    description='A Python implementation of flux balance analysis',
    tests_require = ['nose'],
    cmdclass={'test': NoseTestCommand},
    install_requires=['beautifulsoup4', 'glpk>=0.3.1', 'python-libsbml-experimental'],
    include_package_data=True,
    dependency_links = ['https://github.com/bradfordboyle/pyglpk/tarball/master#egg=glpk-0.3.1'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Unix',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.0',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
)

