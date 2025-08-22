import numpy as np
from scipy.special import sph_harm

def wavefunction(r, θ, φ, k, f_theta=None):
    """
    Quantum wavefunction for scattering
    ψ(r, θ, φ) = e^{ikz} + f(θ) e^{ikr}/r
    
    Parameters:
    r : float - radial distance
    θ : float - polar angle in radians
    φ : float - azimuthal angle in radians
    k : float - wave number
    f_theta : callable - scattering amplitude function f(θ)
    
    Returns:
    complex - wavefunction value
    """
    # Incident wave (e^{ikz} = e^{ikr cosθ})
    incident = np.exp(1j * k * r * np.cos(θ))
    
    # Scattered wave
    if f_theta is not None:
        scattered = f_theta(θ) * np.exp(1j * k * r) / r
    else:
        scattered = 0
    
    return incident + scattered

def scattering_amplitude(theta, k, delta_l=None, potential=None):
    """
    Compute scattering amplitude f(θ)
    
    Options:
    1. From phase shifts: f(θ) = (1/k) Σ (2l+1) e^{iδ_l} sinδ_l P_l(cosθ)
    2. From Born approximation (if potential provided)
    """
    from scipy.special import legendre
    
    if delta_l is not None:
        # Partial wave expansion
        f = 0
        for l, δ in enumerate(delta_l):
            f += (2*l + 1) * np.exp(1j*δ) * np.sin(δ) * legendre(l)(np.cos(theta))
        return f / k
    
    elif potential is not None:
        # Born approximation implementation
        from .born_approximation import first_born_approximation
        return first_born_approximation(theta, k, potential)
    
    else:
        raise ValueError("Must provide either delta_l or potential")