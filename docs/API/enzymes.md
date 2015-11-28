# The [enzymes](../PyFBA/metabolism/enzyme.py) class
 
The [enzymes](../PyFBA/metabolism/enzyme.py) class represents an enzyme, which is the union point of reactions and 
roles.

Each enzyme has a 

* `name`
* a `set` of `roles` associated with it
* a `dict` of `pegs` that are connected to it. The key is the peg, the value is the functional role
* a `dict` of `roles_with_pegs` that are connected. The key is the functional role and the value is a list of pegs with that role
* a `set` of `reactions` that the enzyme connects to

The Enzyme class has the following methods:

