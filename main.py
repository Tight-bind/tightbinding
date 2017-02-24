import numpy as np
from scipy.optimize import bisect
from functools import partial
from get_input_files import TightBindingParameters
from HamiltonianMatrix import HamiltonianMatrix
"""
Author: Jack Baker
Date: 27/10/16

Main script for testiing code.
"""
#for i in range(10):
#    silicon = HamiltonianMatrix(kpoint=np.array([0, 0, 0],dtype=np.float64))
#    print("done")
#silicon = HamiltonianMatrix(kpoint=np.array([0, 0, 0,], dtype=np.float64))
#eigenenergies = silicon.calc_bands(fromfile=True)
#import matplotlib.pyplot as plt
#for i in range(32):
#    plt.plot(eigenenergies[i])
#
#plt.savefig("si_bands.png")
#plt.show()
kb = 1.38064852e-23
class TotalEnergy(TightBindingParameters):

    def __init__(self):
        super(TotalEnergy, self).__init__()

    @staticmethod
    def fermi_dirac(smear_temp,Enk, Ef):
        return np.sum((1 + np.exp((Enk - Ef)/(kb*smear_temp)))**(-1))


    def get_all_eigenvalues(self):
        kpts = ((2*np.pi)/self.simcelldims[0])*self.kpoints
        num_kpt = len(kpts)
        all_energies = np.zeros(num_kpt*self.total_orbitals)
        for ikpt in range(num_kpt):
            H = HamiltonianMatrix(kpts[ikpt]).H
            energy, vector = np.linalg.eigh(H)
            energy = energy.real
            for ien in range(len(energy)):
                print(energy[ien])
                idx = ikpt*32 + ien
                all_energies[idx] = energy[ien]
                print(self.smear_temp)
        obj_func = partial(self.fermi_dirac,
                           smear_temp=self.smear_temp,
                           Enk=all_energies)
        fermi_en = bisect(obj_func, 50, -10)


if __name__ == "__main__":
    en = TotalEnergy()
    en.get_all_eigenvalues()

