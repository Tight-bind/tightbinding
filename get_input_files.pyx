import numpy as np
from cmath import exp
from itertools import product
from Bakerbind_input import *

"""
Author: Jack Baker

Date: 20/02/17

Description: The base class used to retrieve all of the parameters from
             Bakerbind_input and the atomic coordinate file. All input
             parameters are then written to input.log.
             We don't want global variables in HamiltonianMatrix hence why
             Bakerbind_input parameters are set to attributes here.
"""

cdef class TightBindingParameters:

    @staticmethod
    def make_mp_mesh(simcelldims, mpdims):
        # evenly sample unit cube -0.5 < x < 0.5
        spacing = 1/mpdims
        kx_points = np.arange(-0.5 + spacing[0]/2 , 0.5, spacing[0])
        ky_points = np.arange(-0.5 + spacing[1]/2 , 0.5, spacing[1])
        kz_points = np.arange(-0.5 + spacing[2]/2 , 0.5, spacing[2])
        k_points = product(kx_points, ky_points, kz_points)
        # TODO: Work out irreducible BZ points. This is inefficient.
        return np.array([simcelldims*kpoint for kpoint in list(k_points)])

    def __init__(self):
        # Begin with Bakerbind_input
        self.job_name = job_name
        self.coordinate_file = coordinate_file
        self.simcelldims = np.array([sim_box_x, sim_box_y, sim_box_z])
        self.eta_dict = eta_dict
        self.on_site_dict = on_site_dict
        self.sol_method = sol_method
        self.isfractionalcoord = isfractionalcoord
        self.dist_cut_off = dist_cut_off
        # Now coordinate file parameters

        coord_data = np.genfromtxt(coordinate_file,
                                   dtype=['U15', '<f8', '<f8', '<f8'
                                          , '<f8', '<f8'],
                                   names=('Species', 'x_coord',
                                          'y_coord', 'z_coord',
                                          'num_orbitals', 'num_valence'),
                                   skip_header=1)
        self.num_atoms = len(coord_data)
        self.coords = np.array([list(coord)[1:-2] for
                                coord in coord_data],
                                dtype=np.float)
        self.species = np.array([coord[0] for
                                 coord in coord_data],
                                 dtype=np.str)
        self.num_orbitals = np.array([coord[4] for
                                      coord in coord_data],
                                      dtype=np.int) 
        self.valence_charges = np.array([coord[5] for
                                         coord in coord_data],
                                         dtype=np.int)
        self.num_electrons = np.sum(self.valence_charges)
        self.total_orbitals = np.sum(self.num_orbitals)
        if isfractionalcoord:
            self.coords = np.array([coord*self.simcelldims for
                                    coord in self.coords],
                                    dtype=np.float)
        self.mpx = mpx
        self.mpy = mpy
        self.mpz = mpz
        self.kpoints = self.make_mp_mesh(self.simcelldims, np.array([mpx, mpy, mpz]))
        # send everything to input.log
        with open("input.log", "w") as input_log:
            for key, item in self.__dict__.items():
                input_log.write("   =    ".join([str(key), str(item) + "\n"]))

if __name__ == "__main__":
    params = TightBindingParameters()
