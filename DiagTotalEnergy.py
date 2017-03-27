from scipy.optimize import bisect
from functools import partial
from get_input_files import TightBindingParameters
import numpy as np
from HamiltonianMatrix import HamiltonianMatrix
"""
Author: Jack Baker

Date Created: 27/02/17:

Description: calculate the total energy from a diagonalization solution by
bisecting the fermi-dirac distribution under the contraint of the number of
electrons.
"""

class DiagTotalEnergy(TightBindingParameters):

    def __init__(self):
        super(DiagTotalEnergy, self).__init__()
        all_eig_values = self.get_all_eigenvalues()
        obj_func = partial(self.fermi_dirac, all_eig_values)
        fermi_en = bisect(obj_func, all_eig_values.min(),
                                    all_eig_values.max())
        print(fermi_en)


    def fermi_dirac(self, all_eig_values, Ef):
        """
        INPUTS: Ef => Fermi energy
        OUTPUTS: => Non linear equation to bisect for Ef.
        """
        kb_by_ev = 0.86173213995
        return np.sum((1 + np.exp((all_eig_values -\
        Ef)/(kb_by_ev*self.smear_temp)))**(-1)) - self.num_electrons


    def get_all_eigenvalues(self):
        """
        INPUTS: None
        OUTPUTS: All of the energies for each band and K point in an array.
        """
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
        print(all_energies)
        return all_energies

if __name__ == "__main__":
    DiagTotalEnergy()
