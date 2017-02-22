import numpy as np
cimport numpy as np
from get_input_files import TightBindingParameters
from slater_koster import on_site_energy_table, slater_koster_table,\
                          bloch_phase_factor
from math import sqrt

"""
Author: Jack Baker

Date Created: 28/10/16: Modified 21/02/17

Description: Build and solve (by diagonalization) the tight-binding
Hamiltonian matrix using the Slater-Koster tables for the electron hopping
elements.
"""
ctypedef np.complex128_t COMPLEX
ctypedef np.float64_t DOUBLE
ctypedef np.int_t INT



class HamiltonianMatrix(TightBindingParameters):
    """
    A child class of TightBindingParameters used for building and solving
    the tight-binding Hamiltonian matrix.
    """
    def __init__(self, np.ndarray kpoint, int N_images=1):
        super(HamiltonianMatrix, self).__init__()
        cdef int total_orbitals = self.total_orbitals
        cdef np.ndarray[COMPLEX, ndim=2] H = np.zeros([total_orbitals,
                                                           total_orbitals],
                                                           dtype=np.complex)
        cdef np.ndarray[DOUBLE, ndim=1] simcell = self.simcelldims
        cdef int num_atoms = self.num_atoms
        cdef np.ndarray[INT, ndim=1] orbitals = self.num_orbitals
        cdef np.ndarray orb_order = np.array(["ss", "px", "py", "pz"],
                                       dtype=np.str)
        cdef np.ndarray[DOUBLE, ndim=2] coords = self.coords
        cdef double dist_cut_off = self.dist_cut_off
        cdef np.ndarray species = self.species
        cdef np.ndarray[DOUBLE, ndim=1] ri, rj, rij
        #cdspeciesi, orbstri, orbstrj
        cdef double distij
        cdef int iatom, jatom, iorb, jorb, orbidxi, orbidxj
        cdef int imx, imy, imz
        for iatom in range(num_atoms):
          for jatom in range(num_atoms):
            ri, rj = coords[iatom], coords[jatom]
            speciesi = species[iatom]
            for iorb in range(orbitals[iatom]):
              for jorb in range(orbitals[jatom]):
                orbstri, orbstrj = orb_order[iorb], orb_order[jorb]
                orbidxi = np.sum(orbitals[:iatom]) + iorb
                orbidxj = np.sum(orbitals[:jatom]) + jorb
                for imx in range(N_images, N_images+1):
                  for imy in range(N_images, N_images+1):
                    for imz in range(N_images, N_images+1):
                      rij = ri - rj
                      rij[0] = imx*simcell[0] - rij[0]
                      rij[1] = imy*simcell[1] - rij[1]
                      rij[2] = imz*simcell[2] - rij[2]
                      distij = np.linalg.norm(rij)
                      if distij > dist_cut_off:
                        continue
                      elif (iatom==jatom and imx==0 and imy==0 and imz==0):
                        if orbstrj == orbstri:
                          H[orbidxi, orbidxj] = H[orbidxi, orbidxj] +\
                          on_site_energy_table(speciesi,
                                               orbstri)*bloch_phase_factor(kpoint,
                                                                           rij)
                        else:
                          continue
                      else:
                        H[orbidxi, orbidxj] = H[orbidxi, orbidxj] +\
                                              slater_koster_table(
                                               orbstri, orbstrj,
                                               rij[0]/distij,
                                               rij[1]/distij,
                                               rij[2]/distij,
                                               distij)*bloch_phase_factor(kpoint,
                                                                          rij)
        self.H = H

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

    def calc_bands(self, fromfile=False):
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
        e_array = [[] for i in range(self.total_orbitals)]
        if fromfile:
            self.kpts = np.genfromtxt("kpoints.in", delimiter=",")

        else:
            # Default trajectory: W=>G=>X=>W=>L=>G
            # WARNING. HARD-CODED FOR SILICON LATTICE PARAMETER
            # W => G
            self.kpts = [((2*np.pi)/5.431)*np.array([kx, ky, 0]) for
                         kx, ky in zip(np.linspace(0.5, 0, 50),
                         np.linspace(1, 0, 20))]
            # G => X
            self.kpts.extend([((2*np.pi)/5.431)*np.array([0, ky, 0]) for
                         ky in np.linspace(0, 1, 20)])
            # X => W
            self.kpts.extend([((2*np.pi)/5.431)*np.array([kx, 1, 0]) for
                         kx in np.linspace(0, 0.5, 20)])
            # W => L
            self.kpts.extend([((2*np.pi)/5.431)*np.array([0.5, ky, kz]) for
                         ky, kz in zip(np.linspace(1, 0.5, 20),
                         np.linspace(0, 0.5, 50))])
            # L => G
            self.kpts.extend([((2*np.pi)/5.431)*np.array([kx, kx, kx]) for
                         kx in np.linspace(0.5, 0, 20)])
        for kpt in self.kpts:
            self.build_TB_H(kpt)
            energies = self.solve_H()
            for i in range(self.total_orbitals):
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
    bulk_si_H = HamiltonianMatrix(kpoint=np.array([0, 0, 0,]))
    bulk_si_H.show_matrix()
    bulk_si_H = HamiltonianMatrix(kpoint=np.array([0, 0, 0,]))
    bulk_si_H.show_matrix()
    bulk_si_H = HamiltonianMatrix(kpoint=np.array([0, 0, 0,]))
    bulk_si_H.show_matrix()
    #bulk_si_H.matrix_to_csv(file_name="newmethodH.csv")
    #eigenenergies = bulk_si_H.calc_bands()
    #import matplotlib.pyplot as plt
    #for i in range(32):
    #    plt.plot(eigenenergies[i])
    #plt.savefig("bands.png")
