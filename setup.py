import setuptools
from distutils.core import setup
from setuptools.command.test import test as TestCommand

import io

# this is taken from https://www.jeffknupp.com/blog/2013/08/16/open-sourcing-a-python-project-the-right-way/ and
# allows README.txt and CHANGES.txt to automatically be used by PyPI. Create README.txt using http://pandoc.org/:
# pandoc -t plain -o README.txt README.md


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.txt')

def get_version():
    with open("VERSION", 'r') as f:
        v = f.readline().strip()
        return v


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
    version=get_version(),
    packages=setuptools.find_packages(),
    url='http://linsalrob.github.io/PyFBA/',
    license='The MIT License (MIT)',
    author='Rob Edwards',
    author_email='raedwards@gmail.com',
    long_description=long_description,
    platforms='any',
    install_requires=[
        "beautifulsoup4",
        "jupyter",
        "lxml",
        "nose",
        "python-libsbml",
        'importlib_resources; python_version < "3.7"',
        'glpk'
    ],
    test_suite = 'nose.collector',
    description='A Python implementation of flux balance analysis',
    tests_require = ['nose'],
    cmdclass={'test': NoseTestCommand},
    include_package_data=True,
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
    ],
    entry_points={
        "console_scripts": ["pyfba = PyFBA.cmd.entry:run"]
    },
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    zip_safe=True
)

