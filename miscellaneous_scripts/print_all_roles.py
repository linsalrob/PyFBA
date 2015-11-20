import os
import sys

import filters
from parse import model_seed

__author__ = 'Rob Edwards'

"""
Print all roles with reactions
"""

# rxn00867
# Histidine ammonia-lyase (EC 4.3.1.3)



reactions =  filters.roles_to_reactions({'2-ketoaldonate reductase, broad specificity (EC 1.1.1.215) (EC 1.1.1.-)', '4-hydroxybenzoate transporter',
          'ABC transporter permease protein YejB', 'ABC transporter permease protein YejE',
          'Bacterioferritin-associated ferredoxin', 'Cd(II)/Pb(II)-responsive transcriptional regulator',
          'Conidiation-specific protein 10', 'DNA-binding response regulator KdpE',
          'FIG002060: uncharacterized protein YggL', 'FIG00510550: hypothetical protein',
          'FIG00510741: hypothetical protein', 'FIG00638130: hypothetical protein', 'FIG00640809: hypothetical protein',
          'FIG01046673: hypothetical protein', 'Fimbrial protein precursor', 'Fumarate hydratase class I (EC 4.2.1.2)',
          'Gamma-glutamyltranspeptidase (EC 2.3.2.2) @ Glutathione hydrolase (EC 3.4.19.13)',
          'Glycine cleavage system transcriptional activator', 'Hypothetical metal-binding enzyme, YcbL homolog',
          'IS5 transposase', 'NADPH:quinone oxidoreductase 2', 'Probable NreB protein',
          'Protocatechuate 3,4-dioxygenase beta chain (EC 1.13.11.3)', 'Sucrose phosphorylase (EC 2.4.1.7)',
          'Transcriptional regulator of D-allose utilization, RpiR family',
          'Transcriptional repressor protein TrpR # regulates the trp, aroA I (aroH, aroL) and the mtr operons and auto-regulates TrpR.',
          'two component transcriptional regulator, winged helix family',
          'Uncharacterized iron-regulated membrane protein; Iron-uptake factor PiuB',
          'Uncharacterized PLP-dependent aminotransferase YfdZ', 'Ync'})

for r in reactions:
    print("{}\t{}".format(r, reactions[r]))

sys.exit(0)

reactions = filters.roles_to_reactions({'Histidine ammonia-lyase (EC 4.3.1.3)'}, verbose=True)
print("HAL")
for r in reactions:
    print("{}\t{}".format(r, reactions[r]))
print("LAH\n\n")

roles = model_seed.roles()
reactions = filters.roles_to_reactions(set(roles.keys()))

for r in reactions:
    print("{}\t{}".format(r, reactions[r]))