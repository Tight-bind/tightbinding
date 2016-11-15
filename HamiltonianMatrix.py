import numpy as np
from read_coord import AtomicCoordinates
from slater_koster import on_site_energy_table, slater_koster_table

"""
Author: Jack Baker

Date: 28/10/16

Description: Build the tight binding hamiltonian using the
Slater-Koster tables.
"""


class HamiltonianMatrix(AtomicCoordinates):
    """
    A child class of Atomic_Coordinates used for building and solving
    the tight binding Hamiltonian matrix.
    """

    def __init__(self, coord_file, species_dict={}, coords_dict={},
                 orbital_dict={}, atom_id_dict={}):
        AtomicCoordinates.__init__(self, coord_file, species_dict,
                                   coords_dict, orbital_dict, atom_id_dict)
        AtomicCoordinates.generate_dict(self)
        total_orbitals = len(orbital_dict)
        self.H = np.empty([total_orbitals, total_orbitals])


    def calc_matrix_elements(self, show_matrix, dist_cut_off):
        for i, orb_i in self.orbital_dict.items():
            for j, orb_j in self.orbital_dict.items():
                if i == j:
                    element = on_site_energy_table(self.species_dict[
                                                        self.atom_id_dict[i]], orb_i)
                elif orb_i != orb_j and self.atom_id_dict[i] == self.atom_id_dict[j]:
                    element = 0
                else:
                    r_i = self.coords_dict[self.atom_id_dict[i]]
                    r_j = self.coords_dict[self.atom_id_dict[j]]
                    r_ij = r_i - r_j
                    dist_ij = np.sqrt(np.dot(r_ij, r_ij))
                    if dist_ij > dist_cut_off:
                        element = 0
                    else:
                        x_dir_cos = r_ij[0]/dist_ij
                        y_dir_cos = r_ij[1]/dist_ij
                        z_dir_cos = r_ij[2]/dist_ij
                        element = slater_koster_table(orb_i, orb_j, x_dir_cos,
                                                      y_dir_cos, x_dir_cos, dist_ij)
                self.H[i, j] = element
        if show_matrix:
            print("\n****** Hamiltonian Matrix ******\n")
            np.set_printoptions(precision=2)
            print(np.matrix(self.H))


    def solve_H(self):
        eigen_energies, eigen_vectors = np.linalg.eig(self.H)
        print("\n****** Eigenvalues and Eigenvectors ******* \n")
        for idx, eigen_energy in enumerate(eigen_energies):
            print("\nSolution " + str(idx + 1))
            print("Eigenenergy: " + str(eigen_energy) + " eV")
            print("Eigenvector: " + str(eigen_vectors[:,idx]))

if __name__ == "__main__":
    CH4_matrix = HamiltonianMatrix("CH4.coord")
    CH4_matrix.calc_matrix_elements(True, 1.5)
    CH4_matrix.solve_H()

