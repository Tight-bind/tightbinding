import numpy as np
from read_coord import AtomicCoordinates

"""
Author: Jack Baker

Date: 28/10/16

Description: Build the tight binding hamiltonian using the
Slater-Koster tables.
"""

mass_electron = 9.10938356e-31 #kg
reduc_pl_const = 1.054571800e-34


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
    3-tuple: r_ij: the x, y, z coordinates for 
    """
    V = (eta * reduc_pl_const ** 2)/(mass_electron * internuclear_distance ** 2)
    return V


def slater_koster_table(orb_i, orb_j, x_dir_cos, y_dir_cos, z_dir_cos, dist_ij):
    if orb_i == "ss" and orb_j == "ss":
        element = V_coeff(eta_coeff("s", "s", "sigma"), dist_ij)
    elif orb_i == "ss" and orb_j == "px":
        element = x_dir_cos * V_coeff(eta_coeff("s", "p", "sigma"), dist_ij)
    elif orb_i == "px" and orb_j == "px":
        element = x_dir_cos ** 2 * V_coeff(eta_coeff("p", "p", "sigma"), dist_ij) +\
        (1 - x_dir_cos ** 2) * V_coeff(eta_coeff("p", "p", "pi"), dist_ij)
    elif orb_i == "px" and orb_j == "py":
        element = x_dir_cos * y_dir_cos *V_coeff(eta_coeff("p", "p", "sigma"), dist_ij) - \
                   x_dir_cos * y_dir_cos *V_coeff(eta_coeff("p", "p", "pi"), dist_ij)
    elif orb_i == "px" and orb_j == "pz":
        element =  x_dir_cos * z_dir_cos *V_coeff(eta_coeff("p", "p", "sigma"), dist_ij) - \
                      x_dir_cos * z_dir_cos *V_coeff(eta_coeff("p", "p", "pi"), dist_ij)

def on_site_energy_table(species, orbital):
    if species == 'C':
        if orbital == 'px' or orbital == 'py' or orbital == 'pz' or orbital == 'ss':
            energy = 0
    elif species == 'H':
        energy = 1
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


    def calc_matrix_elements(self, dist_cut_off):
        for i, orb_i in self.orbital_dict.items():
            for j, orb_j in self.orbital_dict.items():
                if orb_i == orb_j and self.atom_id_dict[i] == self.atom_id_dict[j]:
                    element = on_site_energy_table(self.species_dict[
                                                        self.atom_id_dict[i]], orb_i)
                elif orb_i == orb_j and self.atom_id_dict[i] != self.atom_id_dict[j]:
                    element = 0
                else:
                    element = 99
                self.H[i, j] = element
             
        print(self.H)

                
if __name__ == "__main__":
    inst = HamiltonianMatrix("CH4.coord")    
    inst.calc_matrix_elements(0.1)

