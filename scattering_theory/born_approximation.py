import numpy as np
from scipy.integrate import nquad
from .constants import ħ

def first_born_approximation(theta, k, V, m=1.0):
    """
    First Born approximation for scattering amplitude
    
    f(θ) = -m/(2πħ²) ∫ e^{-iΔk·r} V(r) d³r
    
    Parameters:
    theta : float - scattering angle (radians)
    k : float - wave number
    V : callable - potential function V(r)
    m : float - particle mass
    
    Returns:
    complex - scattering amplitude
    """
    # Momentum transfer Δk = 2k sin(θ/2)
    k_mag = 2 * k * np.sin(theta/2)
    
    # For spherically symmetric potential
    def integrand(r, theta_r, phi_r):
        # Only r dependence for spherical symmetry
        return V(r) * np.exp(-1j * k_mag * r * np.cos(theta_r)) * r**2 * np.sin(theta_r)
    
    # Perform 3D integral (spherical coordinates)
    r_lim = [0, np.inf]
    theta_lim = [0, np.pi]
    phi_lim = [0, 2*np.pi]
    
    # Use numerical integration (for practical use, provide finite r_max)
    result, _ = nquad(integrand, [r_lim, theta_lim, phi_lim], 
                      opts={'limit': 100, 'epsabs': 1e-6, 'epsrel': 1e-6})
    
    return -m/(2*np.pi*ħ**2) * result