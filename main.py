import numpy as np
from HamiltonianMatrix import HamiltonianMatrix
import itertools
from random import random
"""
Author: Jack Baker
Date: 27/10/16

Main script for testiing code.
"""

CH4 = HamiltonianMatrix("CH4.coord")
CH4.calc_molecule_elements(dist_cut_off=1.4)
CH4.show_matrix()
# CH4 = HamiltonianMatrix("CH4.coord", simcelldim=np.array([100, 100, 100]))
# CH4.calc_periodic_elements(dist_cut_off=1, N_images=1, k_point=np.array([0, 0, 0]))
# CH4.show_matrix()
energies = CH4.solve_H(geteigvecs=True)
total_energy = 2*np.sum(energies[0])
print(total_energy)
eigvecs = energies[1]
# print(eigvecs[0])
density = np.matrix(np.zeros([8, 8]))
for i, j in itertools.product(range(len(eigvecs)), range(len(eigvecs))):
     density[i, j] = random()
print(density)
np.set_printoptions(precision=1)
# print(density)
# print(density.trace())
t_en = 2*(CH4.H*density).trace()
print(t_en)
tol = 0.1
while not(16 - tol < density.trace() < 16 + tol):
    if density.trace() < 16:
        density  = density*density
    elif density.trace() > 16:
        density = 2*density - density*density
    print(density.trace())
