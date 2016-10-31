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


class HamiltonianMatrix(AtomicCoordinates):
    """
    A child class of Atomic_Coordinates used for building and solving
    the tight binding Hamiltonian matrix.
    """
    
    def __init__(self, coord_file, species_dict={}, coords_dict={}):
        AtomicCoordinates.__init__(self, coord_file, species_dict, coords_dict)
        AtomicCoordinates.generate_dict(self)
        orbital_no = 0
        for key, species in self.species_dict.items():
            searchfile = open("species_log.db", "r")
            for line in searchfile:
                if species in line:
                    orbital_no += int(line[2:])
            searchfile.close()
        self.H = np.empty([orbital_no, orbital_no])
        print(self.H)

     
    @staticmethod
    def eta_coeff(l_i, l_j, bond_type):
        """
        string: l_i: The ith atomic orbital angular momentum Q. no (s, p, d etc.)
        string: l_j: The jth atomic orbital angular momentum Q. no (s, p, d etc.)
        string: bond_type: sigma, pi etc.
        """
        eta_string = AO_i + "_" + AO_j + "_" + bond_type
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

    @staticmethod
    def V_coeff(eta, internuclear_distance):
        """
        float: eta: Returned from the function eta_coeff (dict look
               up from solid state table.)
        3-tuple: r_ij: the x, y, z coordinates for 
        """
        V = (eta * reduc_pl_const ** 2)/(mass_electron * internuclear_distance ** 2)
        return V

    @staticmethod
    def slater_koster_table(orb_i, orb_j, x_dir_cos, y_dir_cos, z_dir_cos, dist_ij):
        if orb_i == "s" and orb_j == "s":
            element = V_coeff(eta_coeff("s", "s", "sigma"), dist_ij)
        elif orb_i == "s" and orb_j == "x":
            element = x_dir_cos * V_Coeff(eta_coeff("s", "p", "sigma"), dist_ij)
        elif orb_i == "x" and orb_j == "x":
            element = x_dir_cos ** 2 * V_coeff(eta_coeff("p", "p", "sigma"), dist_ij) +\
            (1 - x_dir_cos ** 2) * V_coeff(eta_coeff("p", "p", "pi"), dist_ij) 
        elif orb_i == "x" and orb_j == "y":
            element = x_dir_cos * y_dir_cos *V_coeff(eta_coeff("p", "p", "sigma"), dist_ij) - \
                      x_dir_cos * y_dir_cos *V_coeff(eta_coeff("p", "p", "pi"), dist_ij)
        elif orb_i == "x" and orb_j == "z":
            element =  x_dir_cos * z_dir_cos *V_coeff(eta_coeff("p", "p", "sigma"), dist_ij) - \
                      x_dir_cos * z_dir_cos *V_coeff(eta_coeff("p", "p", "pi"), dist_ij)
        return element


    def calc_matrix_elements(self, dist_cut_off):
        """
        Calcuate sthe vector r_ij between each atom and computes the magnitude,
        dist_ij.
        """
        loop_counter = 0
        for key_i, coord_i in list(self.coords_dict.items())[loop_counter:]:
            for key_j, coord_j in list(self.coords_dict.items())[loop_counter:]:
                r_ij = coord_i - coord_j
                dist_ij = np.sqrt(np.dot(r_ij, r_ij))
                x_dir_cos = r_ij[0]/dist_ij
                y_dir_cos = r_ij[1]/dist_ij
                z_dir_cos = r_ij[2]/dist_ij

                print(str(key_i) + " to " +  str(key_j) + ": "
                      + str(r_ij) + " " + str(dist_ij))
            loop_counter += 1

if __name__ == "__main__":
    inst = HamiltonianMatrix("CH4.coord")    
    inst.calc_matrix_elements(0.1)

