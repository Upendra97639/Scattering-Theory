import numpy as np
from scattering_theory.quantum import wavefunction
from scattering_theory.constants import ħ

def test_free_particle_wavefunction():
    """Test wavefunction for free particle (no scattering)"""
    r, θ, φ = 10.0, np.pi/4, 0.0
    k = 1.0
    psi = wavefunction(r, θ, φ, k)
    
    # For free particle, wavefunction should be e^{ikz}
    expected = np.exp(1j * k * r * np.cos(θ))
    assert np.isclose(psi, expected, atol=1e-6)

def test_scattering_amplitude():
    """Test scattering amplitude calculation"""
    # Implement test for known cases
    pass