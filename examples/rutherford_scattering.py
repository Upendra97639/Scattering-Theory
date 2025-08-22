from scattering_theory import plot_cross_section
import numpy as np

def rutherford_cross_section(theta, Z1, Z2, E):
    """
    Calculate Rutherford differential cross section
    
    dσ/dΩ = (Z1 Z2 e² / (16π ε0 E))² / sin⁴(θ/2)
    """
    # Constants
    e = 1.602e-19  # Electron charge (C)
    ε0 = 8.854e-12  # Vacuum permittivity (F/m)
    k = 1/(4*np.pi*ε0)  # Coulomb constant
    
    # Calculate cross section
    term1 = (Z1 * Z2 * e**2 * k)**2
    term2 = (4 * E)**2
    term3 = np.sin(theta/2)**4
    return term1 / (term2 * term3)

# Parameters
Z1, Z2 = 2, 79  # Alpha particle and Gold nucleus
E = 5e6  # Energy in eV
theta = np.linspace(0.01, np.pi, 100)  # Avoid θ=0

# Calculate cross section
sigma = rutherford_cross_section(theta, Z1, Z2, E)

# Plot
plt = plot_cross_section(theta, sigma, "Rutherford Scattering Cross Section")
plt.savefig('rutherford_scattering.png')
plt.show()