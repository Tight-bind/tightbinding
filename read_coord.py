import numpy as np

"""
Author: Jack Baker

Date: 27/10/16

Description: Load text from coordinate file into numpy array.
Next, convert array into a dictionary with a uid to label the
particular atom, the species and a tuple with the  x, y and z
coordinates (x, y, z)
"""

orbital_order = ["ss", "px", "py", "pz"]

class AtomicCoordinates(object):
    """
    A base class for reading in the atomic coordinates file and and 
    performing some manipulations on the vectors.
    """    
    def __init__(self, coord_file, species_dict={}, coords_dict={},
                 orbital_dict={}, atom_id_dict={}):
        self.raw_atom_data = np.genfromtxt(coord_file,
                                           dtype=['U15', '<f8', '<f8', '<f8'],
                                           names=('Species', 'x_coord',
                                                  'y_coord', 'z_coord',),
                                           skip_header=1
                                           )
        self.coords_dict = coords_dict
        self.orbital_dict = orbital_dict
        self.species_dict = species_dict
        self.atom_id_dict = atom_id_dict

    def generate_dict(self):
        """
        Inserts the atomic coordinates and species into two dictiionaires, 
        each with a UID.
        """
        orbital_no = 0
        for idx, atom_dat in enumerate(self.raw_atom_data):
            self.coords_dict.update({idx: np.array([atom_dat[1],
                                                    atom_dat[2],
                                                    atom_dat[3]
                                                   ]
                                                  )
                                   }
                                  )
            self.species_dict.update({idx: atom_dat[0]})
            searchfile = open("species_log.db", "r")
            for line in searchfile:
                if atom_dat[0] in line:
                    for orbital in orbital_order[:int(line[2:])]:
                        self.orbital_dict.update({orbital_no: orbital})
                        self.atom_id_dict.update({orbital_no: idx})
                        orbital_no += 1
        print(self.atom_id_dict)
        print(self.orbital_dict)

    def show_atomic_data(self):
        """
        Show the read-in atomic data on screen.
        """
        print("\nInput file interpreted as:\n")
        for key in range(len(self.species_dict)):
            print("UID: " + str(key) + ", Species: " + str(self.species_dict[key])
                  + ", xyz coordinates: " + str(self.coords_dict[key]))

if __name__ == "__main__":
    data = AtomicCoordinates("CH4.coord")
    data.generate_dict()
    data.show_atomic_data()

