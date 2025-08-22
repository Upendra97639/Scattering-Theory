import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import spherical_jn, spherical_yn

def radial_equation(r, u, l, E, V, m=1.0, ħ=1.0):
    """Radial Schrödinger equation for partial wave analysis"""
    k = np.sqrt(2*m*E)/ħ
    dudr = u[1]
    d2udr2 = (2*m/ħ**2 * (V(r) - E) + l*(l+1)/r**2) * u[0]
    return [dudr, d2udr2]

def phase_shift(l, k, V, r_max=100, r_match=10, m=1.0, ħ=1.0):
    """
    Calculate phase shift δ_l for a given angular momentum quantum number l
    
    Parameters:
    l : int - angular momentum quantum number
    k : float - wave number
    V : callable - V(r) function
    r_max : float - maximum radius for numerical integration
    r_match : float - matching radius for asymptotic solution
    
    Returns:
    float - phase shift in radians
    """
    E = (ħ*k)**2/(2*m)
    
    # Solve radial equation numerically
    sol = solve_ivp(
        lambda r, u: radial_equation(r, u, l, E, V, m, ħ),
        [0.01, r_max],
        [0.0, 1.0],  # Initial conditions (u(0)=0, u'(0)=1)
        t_eval=np.linspace(0.01, r_max, 1000),
        method='DOP853'
    )
    
    # Extract solution at matching point
    r_vals = sol.t
    u_vals = sol.y[0]
    idx = np.argmin(np.abs(r_vals - r_match))
    u_match = u_vals[idx]
    du_match = sol.y[1][idx]
    
    # Asymptotic forms
    jl_match = spherical_jn(l, k*r_match)
    nl_match = spherical_yn(l, k*r_match)
    d_jl = k * (spherical_jn(l-1, k*r_match) - (l+1)/r_match * jl_match)
    d_nl = k * (spherical_yn(l-1, k*r_match) - (l+1)/r_match * nl_match)
    
    # Matching to asymptotic solution: u(r) ~ A [cosδ j_l(kr) - sinδ n_l(kr)]
    # Compute derivative at matching point
    ratio = du_match / u_match
    numerator = ratio * jl_match - d_jl
    denominator = ratio * nl_match - d_nl
    
    # Compute phase shift
    δ_l = np.arctan(numerator / denominator)
    return δ_l

def total_cross_section(k, delta_l):
    """
    Total cross-section from phase shifts
    σ = (4π/k²) Σ (2l+1) sin²δ_l
    """
    sigma = 0
    for l, δ in enumerate(delta_l):
        sigma += (2*l + 1) * (np.sin(δ) ** 2)
    return (4 * np.pi / k**2) * sigma