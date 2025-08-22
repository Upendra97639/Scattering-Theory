import numpy as np

def optical_theorem(k, sigma_total):
    """
    Optical theorem: Im f(0) = k σ_total / (4π)
    
    Parameters:
    k : float - wave number
    sigma_total : float - total cross section
    
    Returns:
    complex - forward scattering amplitude (imaginary part)
    """
    return 1j * k * sigma_total / (4 * np.pi)