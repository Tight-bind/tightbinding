import numpy as np

"""
Author: Jack Baker

Date: 27/10/16

Description: Load text from coordinate file into numpy array.
Next, convert array into a dictionary with a uid to label the
particular atom, the species and a tuple with the  x, y and z
coordinates (x, y, z)
"""

class Atomic_Coordinates(object):
    """
    A base class for reading in the atomic coordinates file and and 
    performing some manipulations on the vectors.
    """    
    def __init__(self, coord_file, species_dict={}, coords_dict={}):
        self.raw_atom_data = np.genfromtxt(coord_file,
                                           dtype=['U15', '<f8', '<f8', '<f8'],
                                           names=('Species', 'x_coord',
                                                  'y_coord', 'z_coord',),
                                           skip_header=1
                                           )
        self.species_dict = species_dict
        self.coords_dict = coords_dict

    def generate_dict(self):
        """
        Inserts the atomic coordinates and species into two dictiionaires, 
        each with a UID.
        """
        for idx, atom_dat in enumerate(self.raw_atom_data):
            self.coords_dict.update({idx: np.array([atom_dat[1],
                                                    atom_dat[2],
                                                    atom_dat[3]
                                                   ]
                                                  )
                                   }
                                  )
            self.species_dict.update({idx: atom_dat[0]})

    def show_atomic_data(self):
        """
        Show the read-in atomic data on screen.
        """
        print("\nInput file interpreted as:\n")
        for key in range(len(self.species_dict)):
            print("UID: " + str(key) + ", Species: " + str(self.species_dict[key])
                  + ", xyz coordinates: " + str(self.coords_dict[key]))

if __name__ == "__main__":
    data = Atomic_Coordinates("CH4.coord")
    data.generate_dict()
    data.show_atomic_data()

