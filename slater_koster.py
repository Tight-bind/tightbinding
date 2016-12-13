import numpy as np
from cmath import exp

"""
Author: Jack Baker

Date: 15/11/2016

Description: Some free floating functions and a dictinaory of lambda fucntions used in
the implementation of the Slater-Koster tables.
"""

def eta_coeff(l_i, l_j, bond_type):
    """
    string: l_i: The ith atomic orbital angular momentum Q. no (s, p, d etc.)
    string: l_j: The jth atomic orbital angular momentum Q. no (s, p, d etc.)
    string: b
    """
    eta_string = l_i + "_" + l_j + "_" + bond_type
    #eta_dict = {"s_s_sigma": -1.40,
    #            "s_p_sigma": 1.84,
    #            "p_s_sigma": 1.84,
    #            "p_p_sigma": 3.24,
    #            "p_p_pi": -0.81  }
    # Dave's parameters
    eta_dict = {"s_s_sigma": -1.938,
                "s_p_sigma": 1.745,
                "p_s_sigma": 1.745,
                "p_p_sigma": 3.050,
                "p_p_pi": -1.075}
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
    # Daves suggestion for bulk -> distance dependecne should be 1
    #V = 7.619964162248216*eta
    V = 7.619964162248216 * eta * (1/internuclear_distance) ** 2
    return V



orb_dict = {"ssss": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    V_coeff(eta_coeff("s", "s", "sigma"), dist_ij),
            "sspx": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    x_dir_cos * V_coeff(eta_coeff("s", "p", "sigma"), dist_ij),
            "pxss": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    -x_dir_cos * V_coeff(eta_coeff("s", "p", "sigma"), dist_ij),
            "sspy": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    y_dir_cos * V_coeff(eta_coeff("s", "p", "sigma"), dist_ij),
            "pyss": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    -y_dir_cos * V_coeff(eta_coeff("s", "p", "sigma"), dist_ij),
            "sspz": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    z_dir_cos * V_coeff(eta_coeff("s", "p", "sigma"), dist_ij),
            "pzss": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    -z_dir_cos * V_coeff(eta_coeff("s", "p", "sigma"), dist_ij),
            "pxpx": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    x_dir_cos ** 2 * V_coeff(eta_coeff("p", "p", "sigma"), dist_ij) +\
                    (1 - x_dir_cos ** 2) * V_coeff(eta_coeff("p", "p", "pi"), dist_ij),
            "pzpz": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    z_dir_cos ** 2 * V_coeff(eta_coeff("p", "p", "sigma"), dist_ij) +\
                    (1 - z_dir_cos ** 2) * V_coeff(eta_coeff("p", "p", "pi"), dist_ij),
            "pypy": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    y_dir_cos ** 2 * V_coeff(eta_coeff("p", "p", "sigma"), dist_ij) +\
                    (1 - y_dir_cos ** 2) * V_coeff(eta_coeff("p", "p", "pi"), dist_ij),
            "pxpy": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    (x_dir_cos * y_dir_cos *V_coeff(eta_coeff("p", "p", "sigma"), dist_ij)) -\
                    (x_dir_cos * y_dir_cos *V_coeff(eta_coeff("p", "p", "pi"), dist_ij)),
            "pypx": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    (x_dir_cos * y_dir_cos *V_coeff(eta_coeff("p", "p", "sigma"), dist_ij)) -\
                    (x_dir_cos * y_dir_cos *V_coeff(eta_coeff("p", "p", "pi"), dist_ij)),
            "pypz": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    (z_dir_cos * y_dir_cos *V_coeff(eta_coeff("p", "p", "sigma"), dist_ij)) -\
                    (z_dir_cos * y_dir_cos *V_coeff(eta_coeff("p", "p", "pi"), dist_ij)),
            "pzpy": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    (z_dir_cos * y_dir_cos *V_coeff(eta_coeff("p", "p", "sigma"), dist_ij)) -\
                    (z_dir_cos * y_dir_cos *V_coeff(eta_coeff("p", "p", "pi"), dist_ij)),
            "pzpx": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    (x_dir_cos * z_dir_cos *V_coeff(eta_coeff("p", "p", "sigma"), dist_ij)) - \
                    (x_dir_cos * z_dir_cos *V_coeff(eta_coeff("p", "p", "pi"), dist_ij)),
            "pxpz": lambda x_dir_cos, y_dir_cos, z_dir_cos, dist_ij:
                    (x_dir_cos * z_dir_cos *V_coeff(eta_coeff("p", "p", "sigma"), dist_ij)) - \
                    (x_dir_cos * z_dir_cos *V_coeff(eta_coeff("p", "p", "pi"), dist_ij))}

def slater_koster_table(orb_i, orb_j, x_dir_cos, y_dir_cos, z_dir_cos, dist_ij):
    """
    Calculates the interatomic matrix elements using the Slater-Koster table.
    string: orb_i, orb_j: ss, px, py, pz...
    remaining arguments are float64.
    """
    orb_string = orb_i + orb_j
    if orb_string not in orb_dict:
        print("Orbital overlap %s not recognised" %(orb_string))
    else:
        element = orb_dict.get(orb_string)(x_dir_cos, y_dir_cos, z_dir_cos, dist_ij)
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
    elif species == 'Si':
        if orbital == 'px' or orbital == 'py' or orbital == 'pz':
            #energy = -6.52
            energy = -5.75 # Daves
        elif orbital == "ss":
            #energy = -13.55
            energy = -12.2 # Daves
    return energy



def bloch_phase_factor(k_point, r_ij):
    """
    The Bloch plane wave phase factor.
    """
    return exp(1J*(k_point[0]*r_ij[0] + k_point[1]*r_ij[1] + k_point[2]*r_ij[2]))
