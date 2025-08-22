import numpy as np
import matplotlib.pyplot as plt
from scattering_theory import first_born_approximation, plot_cross_section
from scipy.integrate import quad

def yukawa_potential(r, alpha=0.5, V0=1.0):
    """Yukawa potential (screened Coulomb): V(r) = V0 * exp(-alpha*r)/r"""
    return V0 * np.exp(-alpha*r) / r

def gaussian_potential(r, sigma=1.0, V0=1.0):
    """Gaussian potential: V(r) = V0 * exp(-r²/(2σ²))"""
    return V0 * np.exp(-r**2/(2*sigma**2))

def analytic_yukawa(theta, k, alpha, V0, m=1.0, ħ=1.0):
    """
    Analytic solution for Yukawa potential in Born approximation
    f(θ) = -2mV0/(ħ²α) * 1/(κ² + α²) where κ = 2k sin(θ/2)
    """
    κ = 2 * k * np.sin(theta/2)
    return -2 * m * V0 / (ħ**2 * alpha * (κ**2 + alpha**2))

def analytic_gaussian(theta, k, sigma, V0, m=1.0, ħ=1.0):
    """
    Analytic solution for Gaussian potential in Born approximation
    f(θ) = -m/(√(2π) ħ²) * (2π)^{3/2} V0 σ³ exp(-κ²σ²/2)
    where κ = 2k sin(θ/2)
    """
    κ = 2 * k * np.sin(theta/2)
    prefactor = -m * V0 * sigma**3 * (2*np.pi)**(3/2) / (np.sqrt(2*np.pi) * ħ**2)
    return prefactor * np.exp(-(κ*sigma)**2 / 2)

def run_demo():
    # Physical parameters
    k = 2.0  # wave number
    V0 = 1.0  # potential strength
    alpha = 0.5  # Yukawa screening parameter
    sigma = 1.0  # Gaussian width
    
    # Angular range (0 to 180 degrees, avoid θ=0)
    theta = np.linspace(0.01, np.pi, 100)
    
    # Initialize arrays for results
    f_yukawa = np.zeros_like(theta, dtype=complex)
    f_gaussian = np.zeros_like(theta, dtype=complex)
    f_yukawa_analytic = np.zeros_like(theta, dtype=complex)
    f_gaussian_analytic = np.zeros_like(theta, dtype=complex)
    
    # Calculate scattering amplitudes
    print("Calculating Born approximation for Yukawa potential...")
    for i, t in enumerate(theta):
        f_yukawa[i] = first_born_approximation(t, k, lambda r: yukawa_potential(r, alpha, V0))
        f_yukawa_analytic[i] = analytic_yukawa(t, k, alpha, V0)
    
    print("Calculating Born approximation for Gaussian potential...")
    for i, t in enumerate(theta):
        f_gaussian[i] = first_born_approximation(t, k, lambda r: gaussian_potential(r, sigma, V0))
        f_gaussian_analytic[i] = analytic_gaussian(t, k, sigma, V0)
    
    # Differential cross section: dσ/dΩ = |f(θ)|²
    sigma_yukawa = np.abs(f_yukawa)**2
    sigma_gaussian = np.abs(f_gaussian)**2
    sigma_yukawa_analytic = np.abs(f_yukawa_analytic)**2
    sigma_gaussian_analytic = np.abs(f_gaussian_analytic)**2
    
    # Plot results
    plt.figure(figsize=(14, 10))
    
    # Yukawa potential results
    plt.subplot(2, 2, 1)
    plt.semilogy(np.degrees(theta), sigma_yukawa, 'b-', linewidth=2, label='Numerical')
    plt.semilogy(np.degrees(theta), sigma_yukawa_analytic, 'r--', linewidth=2, label='Analytic')
    plt.xlabel('Scattering Angle (degrees)')
    plt.ylabel('dσ/dΩ')
    plt.title('Yukawa Potential: Differential Cross Section')
    plt.legend()
    plt.grid(True, which="both", ls="--")
    
    # Gaussian potential results
    plt.subplot(2, 2, 2)
    plt.semilogy(np.degrees(theta), sigma_gaussian, 'b-', linewidth=2, label='Numerical')
    plt.semilogy(np.degrees(theta), sigma_gaussian_analytic, 'r--', linewidth=2, label='Analytic')
    plt.xlabel('Scattering Angle (degrees)')
    plt.ylabel('dσ/dΩ')
    plt.title('Gaussian Potential: Differential Cross Section')
    plt.legend()
    plt.grid(True, which="both", ls="--")
    
    # Polar plot for Yukawa
    plt.subplot(2, 2, 3, projection='polar')
    plt.plot(theta, sigma_yukawa, 'b-', linewidth=2, label='Numerical')
    plt.plot(theta, sigma_yukawa_analytic, 'r--', linewidth=2, label='Analytic')
    plt.title('Yukawa Potential: Polar Plot', pad=20)
    plt.legend(loc='upper right')
    plt.grid(True)
    
    # Polar plot for Gaussian
    plt.subplot(2, 2, 4, projection='polar')
    plt.plot(theta, sigma_gaussian, 'b-', linewidth=2, label='Numerical')
    plt.plot(theta, sigma_gaussian_analytic, 'r--', linewidth=2, label='Analytic')
    plt.title('Gaussian Potential: Polar Plot', pad=20)
    plt.legend(loc='upper right')
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig('born_approximation_demo.png', dpi=300)
    plt.show()
    
    # Print forward scattering values
    print("\nResults for forward scattering (θ=0°):")
    print(f"Yukawa potential - Numerical: |f(0)|² = {sigma_yukawa[0]:.4e}")
    print(f"Yukawa potential - Analytic:  |f(0)|² = {sigma_yukawa_analytic[0]:.4e}")
    print(f"Gaussian potential - Numerical: |f(0)|² = {sigma_gaussian[0]:.4e}")
    print(f"Gaussian potential - Analytic:  |f(0)|² = {sigma_gaussian_analytic[0]:.4e}")
    
    # Print backward scattering values
    print("\nResults for backward scattering (θ=180°):")
    print(f"Yukawa potential - Numerical: |f(180°)|² = {sigma_yukawa[-1]:.4e}")
    print(f"Yukawa potential - Analytic:  |f(180°)|² = {sigma_yukawa_analytic[-1]:.4e}")
    print(f"Gaussian potential - Numerical: |f(180°)|² = {sigma_gaussian[-1]:.4e}")
    print(f"Gaussian potential - Analytic:  |f(180°)|² = {sigma_gaussian_analytic[-1]:.4e}")

if __name__ == "__main__":
    run_demo()