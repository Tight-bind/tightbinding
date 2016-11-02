import numpy as np
from read_coord import AtomicCoordinates

"""
Author: Jack Baker

Date: 28/10/16

Description: Build the tight binding hamiltonian using the
Slater-Koster tables.
"""


def eta_coeff(l_i, l_j, bond_type):
    """
    string: l_i: The ith atomic orbital angular momentum Q. no (s, p, d etc.)
    string: l_j: The jth atomic orbital angular momentum Q. no (s, p, d etc.)
    string: b
    """
    eta_string = l_i + "_" + l_j + "_" + bond_type
    eta_dict = {"s_s_sigma": -1.40,
                "s_p_sigma": 1.84,
                "p_s_sigma": 1.84,
                "p_p_sigma": 3.24,
                "p_p_pi": -0.81  }
    if eta_string not in eta_dict:
        print("Unknown bond orbital combination")
    else:
        eta = eta_dict.get(eta_string)
        return eta


def V_coeff(eta, internuclear_distance):
    """
    float: eta: Returned from the function eta_coeff (dict look
           up from solid state table.)
    float: internuclear_distance: enters function in angstroms, V in eV 
    """
    V = 7.619964162248216 * eta * (1/internuclear_distance) ** 2
    return V


def slater_koster_table(orb_i, orb_j, x_dir_cos, y_dir_cos, z_dir_cos, dist_ij):
    """
    Calculates the interatomic matrix elements using the Slater-Koster table.
    string: orb_i, orb_j: ss, px, py, pz...
    remaining arguments are float64.
    """
    if orb_i == "ss" and orb_j == "ss":
        element = V_coeff(eta_coeff("s", "s", "sigma"), dist_ij)
    elif (orb_i == "ss" and orb_j == "px") or (orb_i == "px" and orb_j == "ss"):
        element = x_dir_cos * V_coeff(eta_coeff("s", "p", "sigma"), dist_ij)
    elif orb_i == "px" and orb_j == "px":
        element = x_dir_cos ** 2 * V_coeff(eta_coeff("p", "p", "sigma"), dist_ij) +\
        (1 - x_dir_cos ** 2) * V_coeff(eta_coeff("p", "p", "pi"), dist_ij)
    elif (orb_i == "px" and orb_j == "py") or (orb_i == "py" and orb_j == "px"):
        element = x_dir_cos * y_dir_cos *V_coeff(eta_coeff("p", "p", "sigma"), dist_ij) - \
                   x_dir_cos * y_dir_cos *V_coeff(eta_coeff("p", "p", "pi"), dist_ij)
    elif (orb_i == "px" and orb_j == "pz") or (orb_i == "pz" and orb_j == "px"):
        element =  x_dir_cos * z_dir_cos *V_coeff(eta_coeff("p", "p", "sigma"), dist_ij) - \
                      x_dir_cos * z_dir_cos *V_coeff(eta_coeff("p", "p", "pi"), dist_ij)
    elif (orb_i == "py" and orb_j == "ss") or (orb_i == "ss" and orb_j == "py"):
        element = 0
    elif (orb_i == "pz" and orb_j == "ss") or (orb_i == "ss" and orb_j == "pz"):
        element = 0
    return element


def on_site_energy_table(species, orbital):
    """
    Some onsite energies for seleted atomic species.
    string: species: elemental symbol for the species
    string: orbital: px, py ss...
    """
    if species == 'C':
        if orbital == 'px' or orbital == 'py' or orbital == 'pz':
            energy = -8.97
        elif orbital == "ss":
            energy = -17.52
    elif species == 'H':
        energy = -13.61
    return energy


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
                    if dist_ij < dist_cut_off:
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
            print(self.H)

    
    def solve_H(self):
        eigen_energies, eigen_vectors = np.linalg.eig(self.H)
        print("\n****** Eigenvalues and Eigenvectors ******* \n")
        for idx, eigen_energy in enumerate(eigen_energies):
            print("\nSolution " + str(idx + 1))
            print("Eigenenergy: " + str(eigen_energy) + " eV")
            print("Eigenvector: " + str(eigen_vectors[:,idx]))
        
                
if __name__ == "__main__":
    CH4_matrix = HamiltonianMatrix("CH4.coord")    
    CH4_matrix.calc_matrix_elements(True, 1)
    CH4_matrix.solve_H()

