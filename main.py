import numpy as np
#from scipy.optimize import bisect
#from functools import partial
from get_input_files import TightBindingParameters
from HamiltonianMatrix import HamiltonianMatrix
"""
Author: Jack Baker
Date: 27/10/16

Main script for testiing code.
"""
for i in range(100):
    silicon = HamiltonianMatrix(kpoint=np.array([0, 0, 0],dtype=np.float64))
    print("done")
    #print(silicon.H) 
#silicon = HamiltonianMatrix(kpoint=np.array([0, 0, 0,], dtype=np.float64))
#eigenenergies = silicon.calc_bands(fromfile=True)
#import matplotlib.pyplot as plt
#for i in range(32):
#    plt.plot(eigenenergies[i])
#
#plt.savefig("si_bands.png")
#plt.show()
"""
kb = 1.38064852e-23
ev = 1.60218e-19
rat = 0.86173213995
class TotalEnergy(TightBindingParameters):

    def __init__(self):
        super(TotalEnergy, self).__init__()

    
    def fermi_dirac(self, Ef):
        return np.sum((1 + np.exp((self.get_all_eigenvalues() - Ef)/(rat*self.smear_temp)))**(-1)) - self.num_electrons


    def get_all_eigenvalues(self):
        kpts = ((2*np.pi)/self.simcelldims[0])*self.kpoints
        num_kpt = len(kpts)
        all_energies = np.zeros(num_kpt*self.total_orbitals)
        for ikpt in range(num_kpt):
            H = HamiltonianMatrix(kpts[ikpt]).H
            energy, vector = np.linalg.eigh(H)
            energy = energy.real
            for ien in range(len(energy)):
                idx = ikpt*32 + ien
                all_energies[idx] = energy[ien]
        return all_energies

        #fermi_en = bisect(obj_func, 50, -10)
"""

#if __name__ == "__main__":
    #en = TotalEnergy()
    #import matplotlib.pyplot as plt
    #en_space = np.linspace(-50, 50, 50)
    #values = []
    #for eni in en_space:
    #    value = en.fermi_dirac(eni)
    #    values.append(value)
    #plt.plot(en_space, values)
    #plt.show()
    ##print(en.fermi_dirac(10))
    #print(en.fermi_dirac(-50))
    #print(bisect(en.fermi_dirac, 50, -50))

