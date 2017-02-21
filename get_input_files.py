import numpy as np
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

class TightBindingParameters:
    def __init__(self):
        # Begin with Bakerbind_input
        self.job_name = job_name
        self.coordinate_file = coordinate_file
        self.simcelldims = np.array([sim_box_x, sim_box_y, sim_box_z])
        self.eta_dict = eta_dict
        self.on_site_dict = on_site_dict
        self.sol_method = sol_method
        self.isfractionalcoord = isfractionalcoord
        # Now coordinate file parameters

        coord_data = np.genfromtxt(coordinate_file,
                                   dtype=['U15', '<f8', '<f8', '<f8'],
                                   names=('Species', 'x_coord',
                                          'y_coord', 'z_coord'),
                                   skip_header=1)
        print(coord_data)

if __name__ == "__main__":
    TightBindingParameters()
