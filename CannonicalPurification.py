from HamiltonianMatrix import HamiltonianMatrix
import numpy as np

"""
Author: Jack Baker

Date Created: 20/12/16

An implementation of the canonical purification scheme suggested by Palser and
Manopolous as an extenstion of the original work by Mcweeny.
"""

class McweenySolver(HamiltonianMatrix):
    """
    Initialize with a Hamiltonian matrix to perform a total energy calculation
    with linear scaing.
    """
    def __init__(self, coord_file, dist_cut_off,
                 isperiodic=False, N_images=None, k_point=None,
                 isfractionalcoord=False, simcelldim=None):
        HamiltonianMatrix.__init__(self, coord_file, isfractionalcoord,
                                  simcelldim)
        if isperiodic:
            HamiltonianMatrix.calc_periodic_elements(self, dist_cut_off,
                                                      N_images, k_point)
        else:
            HamiltonianMatrix.calc_molecule_elements(self, dist_cut_off)
        num_orbitals = len(self.orbital_dict)
        num_electrons = 32
        mubar = self.H.trace()/num_orbitals
        # get upper and lower bounds on H
        Hmin = min([self.H[i, i] - abs(self.H[i]).sum() + abs(self.H[i, i])
                   for i in range(num_orbitals)])
        Hmax = max([self.H[i, i] + abs(self.H[i]).sum() - abs(self.H[i, i])
                    for i in range(num_orbitals)])
        # The smallest is the Lambda coefficient
        Lambda = min([num_electrons/(Hmax - mubar),
                          (num_orbitals - num_electrons)/(mubar - Hmin)])
        I = np.matrix(np.identity(num_orbitals))
        # Now we can calculate the optimal intial density matrix
        self.density = ((Lambda/num_orbitals)*(mubar*I - self.H)
                        + (num_electrons/num_orbitals)*I)
        np.set_printoptions(precision=2)
        #print(self.density)

    def get_total_energy(self, tol):
        rho = self.density
        H = self.H
        iteration = 0
        residual = 10e8 # initialize with junk to enter while loop
        print(rho)
        print(rho*rho)
        print(rho*rho*rho)
        print((rho*rho - rho*rho*rho).trace())
        print((rho - rho*rho).trace())
        Cn = ((rho*rho - rho*rho*rho).trace()/(rho - rho*rho).trace())[0,0]
        if np.isnan(Cn) == True:
            Cn = 0
        while residual > tol:
            t_en_old = (rho*H).trace()
            # compute unstable fixed point
            Cn = ((rho*rho - rho*rho*rho).trace()/(rho - rho*rho).trace())[0,0]
            if np.isnan(Cn) == True:
                Cn = 0
            if Cn <= 0.5:
                rho = ((1-2*Cn)*rho + (1+Cn)*rho*rho - rho*rho*rho)/(1-Cn)
            else:
                rho = ((1+Cn)*rho*rho - rho*rho*rho)/Cn
            t_en_new = (rho*H).trace()
            iteration += 1
            print(" ".join(["Iteration", str(iteration), str(t_en_new[0,0])]))
            residual = abs(t_en_new - t_en_old)




if __name__ == "__main__":
    solver = McweenySolver("bulkSi.coord", dist_cut_off=2.7, isperiodic=False,
                            N_images=1, isfractionalcoord=True)
    #np.set_printoptions(precision=2)
    #print(solver.H) 
    solver.get_total_energy(tol=1e-12)
    print(solver.solve_H())
    #solver.show_matrix()
    #print(2*np.sum(solver.solve_H()))
    #print(solver.H)
    #eigen_energies, eigen_vectors = np.linalg.eigh(solver.density)
    #print(eigen_energies)