# ModelSEED data

This is a mirror of the [ModelSEEDDatabase](https://github.com/ModelSEED/ModelSEEDDatabase). 
Everyone who uses PyFBA should cite the [ModelSEEDDatabase](https://modelseed.org/), and we
strive to add reminders elsewhere to ensure that you do!

_License:_
The ModelSEEDDatabase is released under the [Creative Commons Attribution 4.0 International 
License](https://creativecommons.org/licenses/by/4.0/) and we have included a copy of the ModelSEED
license here. You can find [the latest 
license](https://github.com/ModelSEED/ModelSEEDDatabase/blob/master/LICENSE) on their Git repo. (As a reminder,
 `PyFBA` is released under the [MIT license](../../../LICENSE))

_Why copy?_ For a long time we required users to download the entire ModelSEED git repo
for inclusion in `PyFBA` but for usability we are migrating to online installs (pip and conda)
so including the files directly eases that transition. Moreover, the ModelSEED database is stable now, and changes
slowly, so we feel more confident in forking it.

_How do I update my data?_ You _should_ be able to copy the latest ModelSEED files directly to these locations
and easily update. If you do that, please make a [pull request](https://github.com/linsalrob/PyFBA/pulls) to share
the update.

_What is here?_ We have included only the data that we parse for the `PyFBA`. All of that is parsed by the 
[parse/model_seed.py](../../parse/model_seed.py) module. 
