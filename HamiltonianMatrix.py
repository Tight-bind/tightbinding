import numpy as np
from read_coord import AtomicCoordinates
from slater_koster import on_site_energy_table, slater_koster_table,\
                          bloch_phase_factor
from math import sqrt

"""
Author: Jack Baker

Date Created: 28/10/16

Description: Build and solve (by diagonalization) the tight-binding
Hamiltonian matrix using the Slater-Koster tables for the electron hopping
elements.
"""


class HamiltonianMatrix(AtomicCoordinates):
    """
    A child class of Atomic_Coordinates used for building and solving
    the tight binding Hamiltonian matrix.
    """
    def __init__(self, coord_file, kpts=[], isfractionalcoord=False,
                 simcelldim=None, species_dict={}, coords_dict={},
                 orbital_dict={}, atom_id_dict={}):
        AtomicCoordinates.__init__(self, coord_file, species_dict,
                                   coords_dict, orbital_dict, atom_id_dict)
        self.simcelldim = simcelldim
        AtomicCoordinates.generate_dict(self, isfractionalcoord, simcelldim)
        total_orbitals = len(orbital_dict)
        self.kpts = kpts
        self.H = np.zeros([total_orbitals, total_orbitals])

    def calc_molecule_elements(self, dist_cut_off):
        """
        Call this method to build the Hamiltonian matrix for an isolated molecule.
        input: dist_cut_off => The distance in which interatomic matrix elements
                               become zero. float64 or int.
        returns: None.
        """
        orbitals = self.orbital_dict.items()
        for i, orb_i in orbitals:
            for j, orb_j in orbitals:
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
        input: dist_cut_off => The distance in which interatomic matrix elements
                               become zero. float64 or int.
        input: N_images => The number of nearest neighbour cell images to
                          include. Usually, set this to 1 as tight-binding will
                          rarely allow overlaps with far away orbitals. int.
        input: k_point => an array of form [kx, ky, kz]. np.array.
        returns: None.

        """
        simcell = self.simcelldim
        self.H = self.H.astype(complex)
        orbitals = self.orbital_dict.items()
        images = range(-N_images, N_images+1)
        for i, orb_i in orbitals:
            for j, orb_j in orbitals:
                element = 0
                r_i = self.coords_dict[self.atom_id_dict[i]]
                r_j = self.coords_dict[self.atom_id_dict[j]]
                for image_x in images:
                    for image_y in images:
                        for image_z in images:
                            r_ij = r_i - r_j
                            r_ij[0] -= image_x*simcell[0]
                            r_ij[1] -= image_y*simcell[1]
                            r_ij[2] -= image_z*simcell[2]
                            dist_ij = sqrt(r_ij[0]*r_ij[0] +
                                           r_ij[1]*r_ij[1] +
                                           r_ij[2]*r_ij[2])
                            if dist_ij > dist_cut_off:
                                continue
                            elif (self.atom_id_dict[i] == self.atom_id_dict[j]
                                  and image_x==0 and image_y==0 and image_z==0):
                                if orb_i == orb_j:
                                    element += on_site_energy_table(
                                               self.species_dict[
                                               self.atom_id_dict[i]],
                                               orb_i
                                               )*bloch_phase_factor(k_point, r_ij)
                                else:
                                    continue
                            else:
                                element += slater_koster_table(
                                                               orb_i, orb_j,
                                                               r_ij[0]/dist_ij,
                                                               r_ij[1]/dist_ij,
                                                               r_ij[2]/dist_ij,
                                                               dist_ij
                                )*bloch_phase_factor(k_point, r_ij)
                self.H[i, j] = element

    def show_matrix(self):
        """
        Displays the current state of the matrix instant in the terminal window
        input: None.
        returns: None.
        """
        print("\n****** Hamiltonian Matrix ******\n")
        np.set_printoptions(precision=2)
        print(self.H)

    def solve_H(self, show_solns=False, geteigvecs=False):
        """
        Solve the Hamiltonian matrix by diagonalization
        input: show_solns => optionally show the eigenvalues and eigenvectors.
                             boolean. Default: False.
        returns: eigen_energies.real => The eigenvalue with the small complex
                                        residual discarded. Float64.
        """
        eigen_energies, eigen_vectors = np.linalg.eigh(self.H)
        print("\nMatrix diagonalized\n")
        if show_solns:
            for idx, eigen_energy in enumerate(eigen_energies):
                print("\nSolution " + str(idx + 1))
                print("Eigenenergy: " + str(eigen_energy) + " eV")
                print("Eigenvector: " + str(eigen_vectors[:,idx]))
        if geteigvecs:
            return eigen_energies.real, eigen_vectors.real
        else:
            return eigen_energies.real

    def calc_bands(self, dist_cut_off, N_images, fromfile=False):
        """
        Calculate the band structure by retrieving eigenvalues given a trajectory
        of points in momentum-space.
        input: dist_cut_off => The distance in which interatomic matrix elements
                               become zero. float64 or int.
        input: N_images => The number of nearest neighbour cell images to
                           include. Usually, set this to 1 as tight-binding will
                           rarely allow overlaps with far away orbitals. int.
        input: fromfile => Optionally allow for K points to be read in from file.
                           Boolean. Default: False.
        returns: e_array => A list of lists containing a the eigenvalues
                            belonging to the band index n. List.
        """
        e_array = [[] for i in range(len(self.orbital_dict))]
        if fromfile:
            self.kpts = np.genfromtxt("kpoints.in", delimiter=",")

        else:
            # Default trajectory: W=>G=>X=>W=>L=>G
            # WARNING. HARD-CODED FOR SILICON LATTICE PARAMETER
            # W => G
            self.kpts = [((2*np.pi)/5.431)*np.array([kx, ky, 0]) for
                         kx, ky in zip(np.linspace(0.5, 0, 50),
                         np.linspace(1, 0, 50))]
            # G => X
            self.kpts.extend([((2*np.pi)/5.431)*np.array([0, ky, 0]) for
                         ky in np.linspace(0, 1, 50)])
            # X => W
            self.kpts.extend([((2*np.pi)/5.431)*np.array([kx, 1, 0]) for
                         kx in np.linspace(0, 0.5, 50)])
            # W => L
            self.kpts.extend([((2*np.pi)/5.431)*np.array([0.5, ky, kz]) for
                         ky, kz in zip(np.linspace(1, 0.5, 50),
                         np.linspace(0, 0.5, 50))])
            # L => G
            self.kpts.extend([((2*np.pi)/5.431)*np.array([kx, kx, kx]) for
                         kx in np.linspace(0.5, 0, 50)])
        calc_periodic_elements = self.calc_periodic_elements
        solve_H = self.solve_H
        for kpt in self.kpts:
            calc_periodic_elements(dist_cut_off, N_images, kpt)
            energies = solve_H()
            for i in range(len(self.orbital_dict)):
                e_array[i].append(energies[i])
        return e_array


    def matrix_to_csv(self, file_name):
        """
        A debug method to send the Hamiltonian matrix to file.
        input: None.
        returns: None.
        """
        np.savetxt(file_name, self.H,
                   fmt='%.3f', delimiter=",")

    def check_symmetry(self):
        """
        A debug method to check if the Hamiltonian matrix is symmetric by
        examining the difference in the matrix and it's conjugate-transpose. If
        the matrix is Hermitian (It should be), the difference should be null.
        input: None.
        returns: None.
        """
        conjtrans = self.H.conj().T
        self.H -= conjtrans
        for element in self.H.flat:
            if abs(element.real) > 1e-12 or abs(element.imag) > 1e-12:
                print("Non-symmetric matrix")
        self.H += conjtrans

if __name__ == "__main__":
    # Example bulk silicon band structure calculation
    bulk_si = HamiltonianMatrix("bulkSi.coord", isfractionalcoord=True,
                                simcelldim=np.array([5.431, 5.431, 5.431]))
    eigenenergies = bulk_si.calc_bands(dist_cut_off=2.7, N_images=1, fromfile=True)
    import matplotlib.pyplot as plt
    for i in range(32):
        plt.plot(eigenenergies[i])
    plt.savefig("bands.png")
