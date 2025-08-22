import numpy as np

def wave_to_energy(wavelength):
    """
    Convert wavelength to energy (for photons)
    E = hc / λ
    """
    h = 4.135667662e-15  # Planck's constant in eV·s
    c = 3e8  # Speed of light in m/s
    return h * c / wavelength

def energy_to_wave(energy):
    """
    Convert energy to wavelength (for photons)
    λ = hc / E
    """
    h = 4.135667662e-15  # Planck's constant in eV·s
    c = 3e8  # Speed of light in m/s
    return h * c / energy