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
    def __init__(self, kpt):
        super(McweenySolver, self).__init__(kpt)
        print("Hamiltonian Built")
        total_orbitals = self.total_orbitals
        num_electrons = self.num_electrons
        mubar = self.H.trace()/num_electrons
        # get upper and lower bounds on H
        Hmin = min([self.H[i, i] - abs(self.H[i]).sum() + abs(self.H[i, i])
                   for i in range(total_orbitals)])
        Hmax = max([self.H[i, i] + abs(self.H[i]).sum() - abs(self.H[i, i])
                    for i in range(total_orbitals)])
        # The smallest is the Lambda coefficient
        Lambda = min([num_electrons/(Hmax - mubar),
                          (total_orbitals - num_electrons)/(mubar - Hmin)])
        I = np.matrix(np.identity(total_orbitals))
        # Now we can calculate the optimal intial density matrix
        self.density = ((Lambda/total_orbitals)*(mubar*I - self.H)
                        + (num_electrons/total_orbitals)*I)
        print("Initial density built")
        #np.set_printoptions(precision=2)

    def get_total_energy(self, tol):
        rho = np.matrix(self.density.astype(float))
        H = np.matrix(self.H.astype(float))
        #Ignore complex warning here. Matrix made at Gamma so no imag part
        iteration = 0
        residual = 2*tol # initialize with junk to enter while loop
        print("Beginning Palser/Manopolous loop")
        while residual > tol:
            t_en_old = (rho*H).trace()
            # compute unstable fixed point
            Cn = ((rho*rho - rho*rho*rho).trace()/(rho - rho*rho).trace())[0,0]
            if Cn <= 0.5:
                rho = ((1-2*Cn)*rho + (1+Cn)*rho*rho - rho*rho*rho)/(1-Cn)
            else:
                rho = ((1+Cn)*rho*rho - rho*rho*rho)/Cn
            t_en_new = (rho*H).trace()
            iteration += 1
            print(" ".join(["Iteration", str(iteration), str(t_en_new[0,0])]))
            residual = abs(t_en_new - t_en_old)




if __name__ == "__main__":
    solver = McweenySolver(np.array([0, 0, 0], dtype=np.float64))
    #np.set_printoptions(precision=2)
    #print(solver.H) 
    solver.get_total_energy(tol=1e-4)
    #print(solver.solve_H())
    #solver.show_matrix()
    #print(2*np.sum(solver.solve_H()))
    #print(solver.H)
    #eigen_energies, eigen_vectors = np.linalg.eigh(solver.density)
    #print(eigen_energies)
