import numpy as np
import pytest
from scattering_theory.partial_waves import phase_shift, total_cross_section
from scipy.special import spherical_jn, spherical_yn

def test_phase_shift_hard_sphere():
    """
    Test phase shift calculation for hard sphere potential
    Known analytic solution: tanδ_l = j_l(ka)/n_l(ka)
    """
    # Hard sphere potential (infinite wall at r=R)
    def hard_sphere(r, R=1.0):
        return np.inf if r < R else 0.0
    
    # Parameters
    k = 2.0  # wave number
    R = 1.0  # sphere radius
    l_values = [0, 1, 2]  # angular momentum quantum numbers
    
    for l in l_values:
        # Numerical solution
        δ_num = phase_shift(l, k, hard_sphere, r_max=10, r_match=5)
        
        # Analytic solution
        jl = spherical_jn(l, k*R)
        nl = spherical_yn(l, k*R)
        δ_analytic = np.arctan(jl/nl)
        
        # Compare within 1% tolerance
        assert np.isclose(δ_num, δ_analytic, rtol=0.01)

def test_phase_shift_square_well():
    """
    Test phase shift calculation for attractive square well potential
    Known analytic solution for s-wave (l=0)
    """
    # Square well potential parameters
    V0 = 10.0  # well depth
    R = 1.0    # well radius
    
    def square_well(r):
        return -V0 if r < R else 0.0
    
    # Parameters for s-wave
    k = 1.5  # wave number
    l = 0
    
    # Numerical solution
    δ_num = phase_shift(l, k, square_well, r_max=10, r_match=5)
    
    # Analytic solution for s-wave
    k_prime = np.sqrt(k**2 + 2*V0)  # wave number inside well
    analytic = np.arctan((k*np.tan(k_prime*R) - k_prime*np.tan(k*R)) / 
                         (k_prime + k*np.tan(k*R)*np.tan(k_prime*R)))
    
    # Compare within 2% tolerance
    assert np.isclose(δ_num, analytic, rtol=0.02)

def test_total_cross_section():
    """
    Test total cross-section calculation using known phase shifts
    """
    # Create dummy phase shifts
    k = 1.0
    delta_l = [0.1, 0.2, 0.3, 0.4]  # phase shifts for l=0,1,2,3
    
    # Calculate total cross-section
    σ_total = total_cross_section(k, delta_l)
    
    # Manual calculation
    manual_total = 0.0
    for l, δ in enumerate(delta_l):
        manual_total += (2*l + 1) * (np.sin(δ) ** 2)
    manual_total *= 4 * np.pi / k**2
    
    # Compare
    assert np.isclose(σ_total, manual_total)

def test_total_cross_section_low_energy():
    """
    Test total cross-section at low energy where only s-wave contributes
    """
    # For low energy, only l=0 contributes significantly
    k = 0.1  # low wave number
    δ0 = 0.5  # s-wave phase shift
    
    # Calculate total cross-section
    σ_total = total_cross_section(k, [δ0])
    
    # Analytic expression for s-wave only
    analytic = 4 * np.pi / k**2 * np.sin(δ0)**2
    
    # Compare
    assert np.isclose(σ_total, analytic)

def test_phase_shift_high_angular_momentum():
    """
    Test that phase shifts approach zero for high angular momentum
    (Centrifugal barrier should prevent interaction)
    """
    # Simple potential
    def potential(r):
        return 1/r  # Coulomb-like
    
    k = 2.0
    high_l = 10  # high angular momentum
    
    # Calculate phase shift
    δ = phase_shift(high_l, k, potential, r_max=10, r_match=5)
    
    # Should be near zero
    assert np.abs(δ) < 0.01

def test_phase_shift_negative_phase():
    """
    Test that attractive potentials produce negative phase shifts
    """
    # Attractive potential
    def attractive_potential(r):
        return -10.0 if r < 2.0 else 0.0
    
    k = 1.0
    l = 0
    
    # Calculate phase shift
    δ = phase_shift(l, k, attractive_potential, r_max=10, r_match=5)
    
    # Should be negative
    assert δ < 0

def test_phase_shift_positive_phase():
    """
    Test that repulsive potentials produce positive phase shifts
    """
    # Repulsive potential
    def repulsive_potential(r):
        return 10.0 if r < 2.0 else 0.0
    
    k = 1.0
    l = 0
    
    # Calculate phase shift
    δ = phase_shift(l, k, repulsive_potential, r_max=10, r_match=5)
    
    # Should be positive
    assert δ > 0

def test_phase_shift_zero_potential():
    """
    Test that with zero potential, phase shift is zero
    """
    def zero_potential(r):
        return 0.0
    
    k = 1.0
    l = 0
    
    # Calculate phase shift
    δ = phase_shift(l, k, zero_potential, r_max=10, r_match=5)
    
    # Should be zero
    assert np.isclose(δ, 0, atol=1e-4)

def test_total_cross_section_optical_theorem():
    """
    Test that total cross-section satisfies the optical theorem
    Im(f(0)) = k σ_total / (4π)
    """
    from scattering_theory.quantum import scattering_amplitude
    
    # Create phase shifts for a potential
    k = 1.5
    delta_l = [0.2, 0.1, 0.05]  # l=0,1,2
    
    # Calculate total cross-section
    σ_total = total_cross_section(k, delta_l)
    
    # Calculate forward scattering amplitude
    f0 = scattering_amplitude(0, k, delta_l)
    
    # Optical theorem: Im(f(0)) = k σ_total / (4π)
    optical_value = k * σ_total / (4 * np.pi)
    
    # Compare imaginary parts
    assert np.isclose(f0.imag, optical_value, rtol=0.01)

def test_phase_shift_convergence():
    """
    Test that phase shifts converge as r_max increases
    """
    # Square well potential
    def square_well(r):
        return -5.0 if r < 2.0 else 0.0
    
    k = 1.0
    l = 0
    
    # Calculate with increasing r_max
    δ1 = phase_shift(l, k, square_well, r_max=5, r_match=3)
    δ2 = phase_shift(l, k, square_well, r_max=10, r_match=5)
    δ3 = phase_shift(l, k, square_well, r_max=20, r_match=10)
    
    # Differences should decrease
    assert np.abs(δ2 - δ1) > np.abs(δ3 - δ2)

if __name__ == "__main__":
    pytest.main([__file__])