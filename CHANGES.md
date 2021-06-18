# Version 2.5

There have been lots of changes, and I am deprecating the change.md file since we use GitHub. Check out the changes 
there!

# Version 1.2

### Crucial edit to setup.py file
Adding PyFBA.model package to list.

# Version 1.1

### Updated all the code to work with Python 3.
However, this breaks compatibility with Python 2 because of the "errors" keyword in open.

The main changes include relative imports, a new class in the tests/ package to allow deep assertions,
and using this type of construct with the open command to ignore unicode errors in ModelData:

	with open(ssfile, 'r', errors='replace') as sin:

### Added function to obtain Model SEED complexes from a set of roles

* PyFBA/filters/roles_and_complexes.py

### Added gap-fill step to add reactions based on existent EC numbers

* PyFBA/gapfill/ecnumbers.py

### Optimized and fixed bugs during gap-fill reaction minimization process

* PyFBA/gapfill/reaction_minimization.py

### New Model class objects

Model object contains functionality to make common processes easier to facilitate.

* Generate model from RAST annotations
* Save model to hard disk location
* Load model from hard disk location
* Run FBA
* Run FBA and obtain flux values
* Run FBA and obtain flux values according to SEED subsystems
* Run gap-fill
* Remember which reactions were gap-filled reactions and which media

### New iPython Notebooks

* Find a metabolite.ipynb
* PATRIC to FBA.ipynb
* Gap-fill_a_model.ipynb
* Saving_and_loading_a_model.ipynb

# Version 1.0

## Added a CITATION.md

Please cite us if you use PyFBA

### Added the iPython Notebook directory
Example code showing:

importing an SBML file and running FBA
* iPythonNotebooks/Using_an_SBML_model.ipynb

how to build your model from functional roles and how to gap-fill that model on a media
* iPythonNotebook/From_functional_roles_to_gap-filling.ipynb

### Added minimization by accuracy
Removed PyFBA/gapfill/minimize_additional_reactions.py and created reaction_minimization.py

### Example gap-fill pipeline
* example_code/gapfill_from_reactions_multiple_conditions.py

# Version 0.951

### Removed the author tag
It is superflous, that is what versioning is for.

# Version 0.95

Several changes to files:

### Updated the subsystem functions and changed the name to a generic name. The date is handled by versioning!

* PyFBA/Biochemistry/SEED/Subsystems/SS_functions.txt
* PyFBA/Biochemistry/SEED/Subsystems/SS_functions_Oct_2015.txt

### Changing the version number

* PyFBA/__init__.py

### Adding reaction flux information that can be pulled after the sm has been solved

* PyFBA/fba/__init__.py
* PyFBA/fba/fluxes.py

### External reactions

Fixed an issue where we were over-incorporating reactions based on things in the media

* PyFBA/fba/bounds.py
* PyFBA/parse/model_seed.py

### Adding verboseness

* PyFBA/fba/run_fba.py
* PyFBA/gapfill/roles.py
* PyFBA/parse/read_media.py

### Reduced complexity

We only test the right half of the reactions if the left half do not grow, and once we have <10 reactions
to search through we just iterate and knock out each reaction sequentially. It is faster than the shuffle 
approach.

* PyFBA/gapfill/minimize_additional_reactions.py

### Changed the way we read roles

SEED has a notion of multifunctional roles, and we added splitting those functions into roles

* PyFBA/gapfill/subsystem.py

### Gap generation code

Testing each of the individual reactions in a model

* PyFBA/gapgeneration/test_reactions.py
* example_code/test_individual_reactions.py

### Shifted the order of gapfilling

This order makes more logical sense!

* example_code/gapfill_from_reactions.py

### Updated the tests

With all these changes, the tests were not right

* PyFBA/tests/reaction_list.txt
* PyFBA/tests/test_fba.py
* PyFBA/tests/test_suggestions.py


# Version 0.9
Refactored separating roles to functions. SEED has a concept of [multifunctional roles](http://www.nmpdr.org/FIG/Html/SEED_functions.html) and this separates out our roles before we search for them. 

# Version 0.8

# Version 0.7


# Version 0.6
* Removed installation dependencies from setup.py because they break installation! You need to install the dependencies manually
* Refactored the code to remove the os.environ dependencies

# Version 0.5
* Initial release
