# Scattering Theory Python Package

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A computational implementation of quantum scattering theory equations for physics research and education.

## Features
- Quantum wavefunction visualization
- Partial wave analysis
- Born approximation
- Coordinate system transformations
- Cross-section calculations

## Installation

### From GitHub
```bash
pip install git+https://github.com/yourusername/scattering-theory.git
```

### Conda Environment
```bash
conda env create -f environment.yml
conda activate scattering-env
```

## Usage
```python
from scattering_theory import scattering_amplitude, plot_cross_section
import numpy as np

# Calculate Rutherford scattering
theta = np.linspace(0.01, np.pi, 100)
k = 5.0  # wave number
Z1, Z2 = 2, 79  # atomic numbers
E = 5e6  # energy in eV

# Calculate cross-section
sigma = (Z1 * Z2 * 1.44e-9)**2 / (4 * E**2 * np.sin(theta/2)**4)

# Plot
plot_cross_section(theta, sigma, "Rutherford Scattering").show()
```

## Examples
See [examples directory](/examples) for complete usage examples.

## Contributing
Contributions are welcome! Please submit a pull request or open an issue.