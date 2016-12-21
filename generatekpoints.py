# A quick script to generate some K points and save them to kpoints.in
import numpy as np

def generate_k_points(a):
    # W => G
    kpts = [((2*np.pi)/a)*np.array([kx, ky, 0]) for
                   kx, ky in zip(np.linspace(0.5, 0, 50),
                   np.linspace(1, 0, 50))]
    # G => X
    kpts.extend([((2*np.pi)/a)*np.array([0, ky, 0]) for
                   ky in np.linspace(0, 1, 50)])
    # X => W
    kpts.extend([((2*np.pi)/a)*np.array([kx, 1, 0]) for
                   kx in np.linspace(0, 0.5, 50)])
    # W => L
    kpts.extend([((2*np.pi)/a)*np.array([0.5, ky, kz]) for
                   ky, kz in zip(np.linspace(1, 0.5, 50),
                   np.linspace(0, 0.5, 50))])
    # L => G
    kpts.extend([((2*np.pi)/a)*np.array([kx, kx, kx]) for
                   kx in np.linspace(0.5, 0, 50)])
    np.savetxt("kpoints.in", kpts,
               fmt='%.3f', delimiter=",")

if __name__ == "__main__":
    generate_k_points(a=5.431)
    points = np.genfromtxt("kpoints.in", delimiter=",")
    print(points)
