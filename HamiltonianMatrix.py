import numpy as np
from read_coord import AtomicCoordinates
from slater_koster import on_site_energy_table, slater_koster_table, image, bloch_phase_factor

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


    def calc_molecule_elements(self, dist_cut_off):
        """
        Call this method to fill the hailtonian matrix with elements relevant to an
        isolated molecule.
        """
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

    def calc_periodic_elements(self, dist_cut_off, simulation_cell_size,
                               N_images, k_point):
        """
        This method calculates the matrix elements for a periodic solid given a
        point in momentum space.
        """
        self.H = self.H.astype(complex)
        for i, orb_i in self.orbital_dict.items():
            for j, orb_j in self.orbital_dict.items():
                total_element = complex(0, 0)
                for image_x in range(-N_images, N_images+1): # 7 = no images + original image
                    for image_y in range(-N_images, N_images+1): # 7 = no images + original image
                        for image_z in range(-N_images, N_images+1): # 7 = no images + original image
                            r_i = self.coords_dict[self.atom_id_dict[i]]
                            r_j = self.coords_dict[self.atom_id_dict[j]]
                            r_ij = r_i - r_j + np.array([image_x*simulation_cell_size,
                                                         image_y*simulation_cell_size,
                                                         image_z*simulation_cell_size])
                            dist_ij = np.sqrt(np.dot(r_ij, r_ij))
                            if dist_ij > dist_cut_off:
                                element = complex(0, 0)
                            elif dist_ij == 0:
                                if orb_i == orb_j:
                                    element = complex(on_site_energy_table(
                                              self.species_dict[self.atom_id_dict[i]], orb_i), 0)
                                elif orb_i != orb_j:
                                    element = complex(0, 0)
                            else:
                                x_dir_cos = r_ij[0]/dist_ij
                                y_dir_cos = r_ij[1]/dist_ij
                                z_dir_cos = r_ij[2]/dist_ij
                                element = complex(slater_koster_table(orb_i, orb_j, x_dir_cos,
                                                          y_dir_cos, x_dir_cos, dist_ij), 0)
                            total_element += element*bloch_phase_factor(k_point, r_ij)
                self.H[i, j] = total_element

    def show_matrix(self):
        print("\n****** Hamiltonian Matrix ******\n")
        np.set_printoptions(precision=2) 
        print(self.H)


    def solve_H(self):
        eigen_energies, eigen_vectors = np.linalg.eig(self.H)
        #print(eigen_energies.real)
        #return eigen_energies.real

        #print("\n****** Eigenvalues and Eigenvectors ******* \n")
        #for idx, eigen_energy in enumerate(eigen_energies):
            #print("\nSolution " + str(idx + 1))
            #print("Eigenenergy: " + str(eigen_energy.real) + " eV")
            #print("Eigenvector: " + str(eigen_vectors[:,idx].real))
        return eigen_energies.real

if __name__ == "__main__":
    CH4_matrix = HamiltonianMatrix("carb.coord")
    #CH4_matrix.calc_molecule_elements(1.5)
    #CH4_matrix.solve_H()
    #CH4_matrix.show_matrix()
    #CH4_matrix.calc_periodic_elements(15, 20, 1, np.array([0, 0, 0]))
    #CH4_matrix.show_matrix()
    #CH4_matrix.solve_H()
    #periodic_CH4 = HamiltonianMatrix("solid_carbon.coord")
    #periodic_CH4.calc_periodic_elements(5, 2, 1, np.array([np.pi, 0, 0]))
    #periodic_CH4.show_matrix()

    energy1 = []
    energy2 = []
    energy3 = []
    energy4 = []
    energy5 = []
    energy6 = []
    energy7 = []
    energy8 = []
    for i in np.linspace(-np.pi/20, np.pi/20, 50):
        CH4_matrix.calc_periodic_elements(15, 20, 1, np.array([i, 0, 0]))
        energies = CH4_matrix.solve_H()
        energy1.append(energies[0])
        energy2.append(energies[1])
        #energy3.append(energies[2])
        #energy4.append(energies[3])
        #energy5.append(energies[4])
        #energy6.append(energies[5])
        #energy7.append(energies[6])
        #energy8.append(energies[7])
    import matplotlib.pyplot as plt
    plt.plot(energy1)
    plt.plot(energy2)
    #plt.plot(energy3)
    #plt.plot(energy4)
    #plt.plot(energy5)
    #plt.plot(energy6)
    #plt.plot(energy7)
    #plt.plot(energy8)
    plt.show()
    
