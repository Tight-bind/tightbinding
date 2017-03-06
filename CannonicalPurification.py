from HamiltonianMatrix import HamiltonianMatrix
import numpy as np
from Bakerbind_input import sparsity

"""
Author: Jack Baker

Date Created: 20/12/16

An implementation of the canonical purification scheme suggested by Palser and
Manopolous as an extenstion of the original work by Mcweeny.
"""

class OrderNSolver(HamiltonianMatrix):
    """
    Initialize with a Hamiltonian matrix to perform a total energy calculation
    with linear scaing.
    """
    def __init__(self, kpt):
        super(OrderNSolver, self).__init__(kpt)
        print("Hamiltonian Built")
        total_orbitals = self.total_orbitals
        num_basis = int(self.total_orbitals/2)
        num_electrons = self.num_electrons
        mubar = self.H.trace()/num_electrons
        # get upper and lower bounds on H
        Hmin = min([self.H[i, i] - abs(self.H[i]).sum() + abs(self.H[i, i])
                   for i in range(num_basis)])
        Hmax = max([self.H[i, i] + abs(self.H[i]).sum() - abs(self.H[i, i])
                    for i in range(num_basis)])
        I = np.matrix(np.identity(num_basis))
        if self.sol_method == "PM":
            # The smallest is the Lambda coefficient
            Lambda = min([num_electrons/(Hmax - mubar),
                              (total_orbitals - num_electrons)/(mubar - Hmin)])
            # Now we can calculate the optimal intial density matrix
            self.density = ((Lambda/total_orbitals)*(mubar*I - self.H)
                            + (num_electrons/total_orbitals)*I)
            print("Initial PM density built")
            #np.set_printoptions(precision=2)
        elif self.sol_method == "SP2":
            self.density = (Hmax*I- self.H)/(Hmax - Hmin)
        elif self.sol_method == "HPCP":
            pass
        elif self.sol_method == "LNV":
            pass
        else:
            print("Unknown solution method")


    def get_PM_energy(self, tol=1e-5):
        rho = np.matrix(self.density.astype(float))
        H = np.matrix(self.H.astype(float))
        #Ignore complex warning here. Matrix made at Gamma so no imag part
        iteration = 0
        residual = 2*tol # initialize with junk to enter while loop
        print("Beginning Palser/Manopolous loop")
        while residual > tol:
            t_en_old = (rho*H).trace()
            rho2 = rho*rho
            rho3 = rho*rho*rho
            # compute unstable fixed point
            Cn = ((rho2 - rho3).trace()/(rho - rho2).trace())[0,0]
            if Cn <= 0.5:
                rho = ((1-2*Cn)*rho + (1+Cn)*rho2 - rho3)/(1-Cn)
            else:
                rho = ((1+Cn)*rho2 - rho3)/Cn
            t_en_new = (rho*H).trace()
            iteration += 1
            print(" ".join(["Iteration", str(iteration), str(t_en_new[0,0])]))
            residual = abs(t_en_new - t_en_old)

    def get_SP2_energy(self, tol=1e-5):
        if sparsity == "sparse":
            from scipy import sparse
            rho = sparse.csc_matrix(self.density.astype(float))
            H = sparse.csc_matrix(self.H.astype(float))
            num_electrons = self.num_electrons
            trace_rho = rho.diagonal().sum()
            residual = 2*tol
            iteration = 0
            while residual > tol:
                t_en_old =(rho*H).diagonal().sum()
                rho_tmp = -rho*rho + rho
                trace_rho_tmp = rho_tmp.diagonal().sum()
                if (abs(2*trace_rho - 2*trace_rho_tmp - num_electrons) >
                    abs(2*trace_rho + 2*trace_rho_tmp - num_electrons)):
                    rho += rho_tmp
                    trace_rho += trace_rho_tmp
                else:
                    rho -= rho_tmp
                    trace_rho -= trace_rho_tmp
                t_en_new =(rho*H).diagonal().sum()
                iteration += 1
                print(" ".join(["Iteration", str(iteration), str(t_en_new)]))
                residual = abs(t_en_old - t_en_new)

        elif sparsity == "dense":
            rho = np.matrix(self.density.astype(float))
            H = np.matrix(self.H.astype(float))
            num_electrons = self.num_electrons
            trace_rho = rho.trace()
            residual = 2*tol
            iteration = 0
            while residual > tol:
                t_en_old =(rho*H).trace()
                rho_tmp = -rho*rho + rho
                trace_rho_tmp = rho_tmp.trace()
                if (abs(2*trace_rho - 2*trace_rho_tmp - num_electrons) >
                    abs(2*trace_rho + 2*trace_rho_tmp - num_electrons)):
                    rho += rho_tmp
                    trace_rho += trace_rho_tmp
                else:
                    rho -= rho_tmp
                    trace_rho -= trace_rho_tmp
                t_en_new =(rho*H).trace()
                iteration += 1
                print(" ".join(["Iteration", str(iteration), str(t_en_new[0,0])]))
                residual = abs(t_en_old - t_en_new)



if __name__ == "__main__":
    solver = OrderNSolver(np.array([0, 0, 0], dtype=np.float64))
    #np.set_printoptions(precision=2)
    #print(solver.H)
    solver.get_SP2_energy()
    #print(solver.solve_H())
    #solver.show_matrix()
    #print(2*np.sum(solver.solve_H()))
    #print(solver.H)
    #eigen_energies, eigen_vectors = np.linalg.eigh(solver.density)
    #print(eigen_energies)
