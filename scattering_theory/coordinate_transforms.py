import numpy as np
from scipy.optimize import root_scalar

def lab_to_cm_angle(theta_lab, m1, m2):
    """
    Convert scattering angle from lab frame to CM frame
    
    tanθ_L = sinθ_C / (cosθ_C + γ)
    γ = m1/m2
    """
    gamma = m1 / m2
    
    # Solve for θ_C using numerical root finding
    def equation(theta_c):
        num = np.sin(theta_c)
        den = np.cos(theta_c) + gamma
        return np.tan(theta_lab) - num/den
    
    # Bracket between 0 and π
    result = root_scalar(equation, bracket=[0, np.pi], method='brentq')
    return result.root

def cross_section_lab_to_cm(sigma_lab, theta_lab, m1, m2):
    """
    Convert cross-section from lab frame to CM frame
    
    σ_L(θ_L) = [ (1 + γ² + 2γ cosθ_C)^{3/2} / |1 + γ cosθ_C| ] σ_C(θ_C)
    """
    gamma = m1 / m2
    theta_cm = lab_to_cm_angle(theta_lab, m1, m2)
    
    numerator = (1 + gamma**2 + 2*gamma*np.cos(theta_cm))**1.5
    denominator = np.abs(1 + gamma*np.cos(theta_cm))
    
    return sigma_lab * numerator / denominator