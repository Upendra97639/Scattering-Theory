Scattering Theory Documentation
Overview
The Scattering Theory Python package provides a comprehensive computational implementation of quantum scattering theory equations for physics research and education. This package enables researchers and students to:

Calculate quantum scattering amplitudes and cross-sections

Visualize wavefunctions and differential cross-sections

Transform between laboratory and center-of-mass reference frames

Implement both partial wave analysis and Born approximation methods

Explore fundamental quantum scattering phenomena

Key Features
Core Functionality
Quantum Wavefunction Calculation: Compute wavefunctions for scattering problems

Partial Wave Analysis: Implement phase shift calculations and partial wave expansions

Born Approximation: First Born approximation for scattering amplitudes

Coordinate Transforms: Convert between lab and CM frames

Cross-Section Calculations: Differential and total cross-sections

Optical Theorem: Implement the fundamental optical theorem

Visualization Tools
Polar plots of differential cross-sections

3D visualizations of quantum wavefunctions

Interactive plotting capabilities

Utilities
Physical constants database

Energy-wavelength conversions

Special function implementations

Installation
From GitHub
bash
pip install git+https://github.com/yourusername/scattering-theory.git
Using Conda
bash
conda env create -f environment.yml
conda activate scattering-env
Quick Start
Basic Usage
python
from scattering_theory import scattering_amplitude, plot_cross_section
import numpy as np

# Calculate scattering amplitude using Born approximation
theta = np.pi/3  # 60 degrees
k = 2.0  # wave number

def coulomb_potential(r):
    return 1/r  # Simple Coulomb potential

f_theta = scattering_amplitude(theta, k, potential=coulomb_potential)

print(f"Scattering amplitude at θ={theta:.2f} rad: {f_theta:.4f}")
Visualizing a Cross-Section
python
import numpy as np
from scattering_theory import plot_cross_section

# Create angular range
theta = np.linspace(0.01, np.pi, 100)

# Rutherford cross-section
Z1, Z2 = 2, 79  # Alpha particle and Gold nucleus
E = 5e6  # Energy in eV
sigma = (Z1 * Z2 * 1.44e-9)**2 / (4 * E**2 * np.sin(theta/2)**4)

# Plot
plt = plot_cross_section(theta, sigma, "Rutherford Scattering")
plt.show()
API Reference
Quantum Module
python
wavefunction(r, θ, φ, k, f_theta=None)
Computes the quantum wavefunction for scattering problems

Parameters:

r: Radial distance

θ: Polar angle in radians

φ: Azimuthal angle in radians

k: Wave number

f_theta: Scattering amplitude function (optional)

Returns: Complex wavefunction value

python
scattering_amplitude(theta, k, delta_l=None, potential=None)
Calculates scattering amplitude using either partial wave analysis or Born approximation

Parameters:

theta: Scattering angle

k: Wave number

delta_l: List of phase shifts (for partial wave analysis)

potential: Potential function V(r) (for Born approximation)

Returns: Complex scattering amplitude

Partial Wave Analysis
python
phase_shift(l, k, V, r_max=100, r_match=10, m=1.0, ħ=1.0)
Calculates phase shift δ_l for given angular momentum quantum number

Parameters:

l: Angular momentum quantum number

k: Wave number

V: Potential function V(r)

r_max: Maximum radius for integration

r_match: Matching radius for asymptotic solution

Returns: Phase shift in radians

python
total_cross_section(k, delta_l)
Calculates total cross-section from phase shifts

Parameters:

k: Wave number

delta_l: List of phase shifts

Returns: Total cross section

Born Approximation
python
first_born_approximation(theta, k, V, m=1.0)
First Born approximation for scattering amplitude

Parameters:

theta: Scattering angle

k: Wave number

V: Potential function V(r)

m: Particle mass

Returns: Complex scattering amplitude

Coordinate Transforms
python
lab_to_cm_angle(theta_lab, m1, m2)
Converts scattering angle from lab frame to CM frame

Parameters:

theta_lab: Scattering angle in lab frame

m1: Mass of incident particle

m2: Mass of target particle

Returns: Scattering angle in CM frame

python
cross_section_lab_to_cm(sigma_lab, theta_lab, m1, m2)
Converts cross-section from lab frame to CM frame

Parameters:

sigma_lab: Differential cross-section in lab frame

theta_lab: Scattering angle in lab frame

m1: Mass of incident particle

m2: Mass of target particle

Returns: Differential cross-section in CM frame

Plotting
python
plot_cross_section(theta, sigma, title="Differential Cross Section")
Creates polar plot of differential cross section

Parameters:

theta: Array of angles

sigma: Array of cross-section values

title: Plot title

Returns: Matplotlib plot object

python
plot_wavefunction_3d(psi_func, r_max=10, k=1.0, n_points=50)
Creates 3D visualization of wavefunction magnitude

Parameters:

psi_func: Wavefunction function

r_max: Maximum radial distance

k: Wave number

n_points: Number of grid points

Returns: Matplotlib 3D plot

Examples
Partial Wave Analysis
python
from scattering_theory import phase_shift, scattering_amplitude, plot_cross_section
import numpy as np

# Define a square well potential
def square_well(r, V0=10, R=1):
    return -V0 if r < R else 0

# Parameters
k = 1.5  # Wave number
max_l = 4  # Maximum angular momentum

# Calculate phase shifts
delta_l = [phase_shift(l, k, square_well) for l in range(max_l+1)]

# Calculate differential cross section
theta = np.linspace(0, np.pi, 100)
sigma = [np.abs(scattering_amplitude(t, k, delta_l))**2 for t in theta]

# Plot
plot_cross_section(theta, sigma, "Square Well Scattering").show()
Born Approximation
python
from scattering_theory import first_born_approximation, plot_cross_section
import numpy as np

# Yukawa potential
def yukawa_potential(r, alpha=0.2):
    return np.exp(-alpha*r)/r

# Calculate scattering amplitudes
theta = np.linspace(0.01, np.pi, 50)
k = 2.0
f_theta = [first_born_approximation(t, k, yukawa_potential) for t in theta]
sigma = np.abs(f_theta)**2

# Plot
plot_cross_section(theta, sigma, "Yukawa Potential Scattering").show()
Contributing
We welcome contributions to the Scattering Theory package! To contribute:

Fork the repository on GitHub

Create a new branch for your feature or bug fix

Implement your changes with appropriate tests

Submit a pull request with a clear description of your changes

Please ensure your code follows PEP 8 style guidelines and includes docstrings for all public functions.

License
This project is licensed under the MIT License - see the LICENSE file for details.

Support
For support, questions, or feature requests, please open an issue on the GitHub repository.