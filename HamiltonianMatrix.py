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

    def __init__(self, coord_file, kpts=[], isfractionalcoord=False, simcelldimensions=None,
                 species_dict={}, coords_dict={}, orbital_dict={},
                 atom_id_dict={}):
        AtomicCoordinates.__init__(self, coord_file, species_dict,
                                   coords_dict, orbital_dict, atom_id_dict)
        self.simcelldimensions = simcelldimensions
        AtomicCoordinates.generate_dict(self, isfractionalcoord, simcelldimensions)
        total_orbitals = len(orbital_dict)
        self.kpts = kpts
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
                                                      y_dir_cos, z_dir_cos, dist_ij)
                self.H[i, j] = element

    def calc_periodic_elements(self, dist_cut_off, N_images, k_point):
        """
        This method calculates the matrix elements for a periodic solid given a
        point in momentum space.
        """
        simcell = np.array([5.431, 5.431, 5.431])
        self.H = self.H.astype(complex)
        for i, orb_i in self.orbital_dict.items():
            for j, orb_j in self.orbital_dict.items():
                total_element = complex(0, 0)
                for image_x in range(-N_images, N_images+1):
                    for image_y in range(-N_images, N_images+1):
                        for image_z in range(-N_images, N_images+1):
                            r_i = self.coords_dict[self.atom_id_dict[i]]
                            r_j = self.coords_dict[self.atom_id_dict[j]]
                            r_ij = r_i - r_j - np.array([image_x*simcell[0],
                                                         image_y*simcell[1],
                                                         image_z*simcell[2]])
                            dist_ij = np.sqrt(np.dot(r_ij, r_ij))
                            if dist_ij > dist_cut_off:
                                element = complex(0, 0)
                            elif orb_i == orb_j and self.atom_id_dict[i] == self.atom_id_dict[j] and image_x==0 and image_y==0 and image_z==0:
                                element = complex(on_site_energy_table(
                                          self.species_dict[self.atom_id_dict[i]], orb_i), 0)
                            elif orb_i != orb_j and self.atom_id_dict[i] == self.atom_id_dict[j] and image_x==0 and image_y==0 and image_z==0:
                                element = complex(0, 0)
                            else:
                                x_dir_cos = r_ij[0]/dist_ij
                                y_dir_cos = r_ij[1]/dist_ij
                                z_dir_cos = r_ij[2]/dist_ij
                                element = complex(slater_koster_table(orb_i, orb_j, x_dir_cos,
                                                          y_dir_cos, z_dir_cos, dist_ij), 0)
                            total_element += element*bloch_phase_factor(k_point, r_ij)
                self.H[i, j] = total_element


    def show_matrix(self):
        print("\n****** Hamiltonian Matrix ******\n")
        np.set_printoptions(precision=2)
        print(self.H)


    def solve_H(self):
        eigen_energies, eigen_vectors = np.linalg.eigh(self.H)
        #print(eigen_energies.real)
        #return eigen_energies.real

        print("\nSolving\n")
        for idx, eigen_energy in enumerate(eigen_energies):
            if abs(eigen_energy.imag) > 1e-12:
                print("Warning, imaginary eigenvalue!")
        #    print("\nSolution " + str(idx + 1))
        #    print("Eigenenergy: " + str(eigen_energy) + " eV")
        #    print("Eigenvector: " + str(eigen_vectors[:,idx]))
        return eigen_energies.real

    def calc_bands(self, dist_cut_off, N_images):
        e_array = [[] for i in range(len(self.orbital_dict))]
        #self.kpts = np.genfromtxt("kpoints.in", skip_header=2)
        self.kpts = [np.pi/5.431*np.array([kx, kx, 0]) for
                     kx in np.linspace(0, 0.5, 50)]
        for kpt in self.kpts:
            self.calc_periodic_elements(dist_cut_off, N_images, kpt)
            self.check_symmetry()
            energies = self.solve_H()
            for i in range(len(self.orbital_dict)):
                e_array[i].append(energies[i])
        return e_array


    def matrix_to_csv(self, file_name):
       np.savetxt(file_name, self.H,
                   fmt='%.3f', delimiter=",")

    def check_symmetry(self):
        conjtrans = self.H.conj().T
        self.H -= conjtrans
        for element in self.H.flat:
            if abs(element.real) > 1e-12 or abs(element.imag) > 1e-12:
                print("Non-symmetric matrix")
        self.H += conjtrans

if __name__ == "__main__":
    bulk_si = HamiltonianMatrix("bulkSi.coord", isfractionalcoord=True,
                                simcelldimensions=np.array([5.431, 5.431, 5.431]))
    #bulk_si.calc_molecule_elements(dist_cut_off=5)
    #bulk_si.calc_periodic_elements(dist_cut_off=2.7, N_images=2, k_point=np.array([0, 0, 0]))
    #bulk_si.solve_H()
    #bulk_si.matrix_to_csv("beforesym.csv")
    #bulk_si.show_matrix()
    #bulk_si.check_symmetry()
    #bulk_si.show_matrix()
    #bulk_si.matrix_to_csv("aftersym.csv")
    #bulk_si.show_matrix()
    #bulk_si.calc_periodic_elements(dist_cut_off=11, N_images=2, k_point=np.array([0, 0, 0]))
    eigenenergies = bulk_si.calc_bands(dist_cut_off=11, N_images=2)
    #bulk_si.calc_molecule_elements(dist_cut_off=11)
    #bulk_si.matrix_to_csv()
    #eigenenergies= bulk_si.calc_bands(dist_cut_off=11, N_images=2)
    import matplotlib.pyplot as plt
    for i in range(32):
        plt.plot(eigenenergies[i])
    plt.show()
