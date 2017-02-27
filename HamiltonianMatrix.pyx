import numpy as np
from HamiltonianMatrix cimport norm
cimport numpy as np
from get_input_files import TightBindingParameters
from slater_koster cimport on_site_energy_table, slater_koster_table,\
                          bloch_phase_factor
from libc.math cimport sqrt
cimport cython

"""
Author: Jack Baker

Date Created: 28/10/16: Modified 21/02/17

Description: Build and solve (by diagonalization) the tight-binding
Hamiltonian matrix using the Slater-Koster tables for the electron hopping
elements.
"""

cdef np.ndarray orb_order = np.array(["ss", "px", "py", "pz"], dtype=np.str)


@cython.boundscheck(False) # turn off bounds-checking
@cython.wraparound(False)  # turn off negative index wrapping
class HamiltonianMatrix(TightBindingParameters):
    """
    A child class of TightBindingParameters used for building and solving
    the tight-binding Hamiltonian matrix.
    """
    def __init__(self, const double [:] kpoint, const int N_images=1):
        super(HamiltonianMatrix, self).__init__()
        cdef int total_orbitals = self.total_orbitals
        cdef complex [:, :] H = np.zeros([total_orbitals,
                                          total_orbitals],
                                          dtype=np.complex)
        cdef double [:] simcell = self.simcelldims
        cdef int num_atoms = self.num_atoms
        cdef long[:] orbitals = self.num_orbitals
        cdef double [:, :] coords = self.coords
        cdef double dist_cut_off = self.dist_cut_off
        species = self.species
        cdef double [:] ri
        cdef double [:] rj
        cdef double [:] rij = np.array([0, 0, 0], dtype=np.float64)
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
                orbidxi = np.sum(orbitals[:iatom]) + iorb # SLOW
                orbidxj = np.sum(orbitals[:jatom]) + jorb # SLOW
                for imx in range(-N_images, N_images+1):
                  for imy in range(-N_images, N_images+1):
                    for imz in range(-N_images, N_images+1):
                      rij[0] = ri[0] - rj[0] - imx*simcell[0]
                      rij[1] = ri[1] - rj[1] - imy*simcell[1]
                      rij[2] = ri[2] - rj[2] - imz*simcell[2]
                      distij = norm(rij)
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
        self.H = np.asarray(H)

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
            self.kpts = np.genfromtxt("kpoints.in", delimiter=",", dtype=np.float64)
            print(self.kpts)

        else:
            print("Not yet supported")
            return 0
        for kpt in self.kpts:
            self.H = self.__class__(kpt).H
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
