Gapfilling
==========

Gap filling is more of an art form than a science!

The goal of gap filling is to add reactions to your metabolic network that complete the network, but not to add so many
reactions that the model will grow under any circumstance at any time. You just want the model to grow in the right
conditions!

PyFBA includes several different gap filling approaches, and the API is designed so that it is easy to add and test your
own gap filling designs.

First, lets take a look at the built in approaches:

.. automodule:: PyFBA.gapfill.essentials
   :members:

.. automodule:: PyFBA.gapfill.maps_to_proteins
   :members:

.. automodule:: PyFBA.gapfill.media
   :members:

.. automodule:: PyFBA.gapfill.orphan_compound
   :members:


.. automodule:: PyFBA.gapfill.probability
   :members:

.. automodule:: PyFBA.gapfill.roles
   :members:

.. automodule:: PyFBA.gapfill.subsystem
   :members:




