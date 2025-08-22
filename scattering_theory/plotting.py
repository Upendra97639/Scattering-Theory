import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_cross_section(theta, sigma, title="Differential Cross Section"):
    """
    Polar plot of differential cross section
    """
    plt.figure(figsize=(8, 8))
    ax = plt.subplot(111, projection='polar')
    ax.plot(theta, sigma, lw=2)
    ax.set_title(title, pad=20)
    ax.grid(True)
    plt.tight_layout()
    return plt

def plot_wavefunction_3d(psi_func, r_max=10, k=1.0, n_points=50):
    """
    3D visualization of wavefunction magnitude
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Create grid in spherical coordinates
    theta = np.linspace(0, 2*np.pi, n_points)
    r = np.linspace(0.1, r_max, n_points)
    R, Theta = np.meshgrid(r, theta)
    
    # Convert to Cartesian coordinates
    X = R * np.cos(Theta)
    Y = R * np.sin(Theta)
    
    # Calculate wavefunction magnitude
    Z = np.zeros_like(X)
    for i in range(n_points):
        for j in range(n_points):
            Z[i, j] = np.abs(psi_func(R[i,j], Theta[i,j], 0, k))
    
    # Plot surface
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', 
                          rstride=1, cstride=1, alpha=0.8)
    
    fig.colorbar(surf, shrink=0.5, aspect=5)
    ax.set_title("Quantum Wavefunction $|\psi(r,\\theta)|$")
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('$|\psi|$')
    plt.tight_layout()
    return plt