from scattering_theory import phase_shift, total_cross_section, scattering_amplitude, plot_cross_section
import numpy as np

# Define a square well potential
def square_well(r, V0=10, R=1):
    return -V0 if r < R else 0

# Parameters
k = 1.5  # Wave number
max_l = 4  # Maximum angular momentum to consider

# Calculate phase shifts
delta_l = []
for l in range(max_l + 1):
    δ = phase_shift(l, k, square_well, r_max=20, r_match=15)
    delta_l.append(δ)
    print(f"Phase shift δ_{l} = {δ:.4f} rad")

# Calculate total cross section
sigma_total = total_cross_section(k, delta_l)
print(f"Total cross section: {sigma_total:.4f}")

# Calculate differential cross section
theta = np.linspace(0, np.pi, 100)
f_theta = [scattering_amplitude(t, k, delta_l) for t in theta]
sigma = np.abs(f_theta)**2

# Plot
plt = plot_cross_section(theta, sigma, "Partial Wave Analysis")
plt.savefig('partial_wave_demo.png')
plt.show()