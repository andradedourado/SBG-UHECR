from scipy.integrate import quad
from scipy.special import zeta
import numpy as np
import scipy.constants as const

m_p = 0.9383e9 # eV
r0 = (const.e**2) / (4 * const.pi * const.epsilon_0 * const.m_e * const.c**2)

# Overflow encountered in exp

# ----------------------------------------------------------------------------------------------------
def photon_density(eps, z, T = 2.7): # Number of photons per unit volume and energy

    return 8 * np.pi / (const.h * const.c)**3 * eps**2 / (np.exp(eps/(const.k * T * (1 + z))) - 1) * (1 + z)**3

# ----------------------------------------------------------------------------------------------------
def phi(kappa):

    c1 = 0.8048
    c2 = 0.1459
    c3 = 1.137 * 1.e-3
    c4 = -3.879 * 1.e-6
    c = [c1, c2, c3, c4]

    d0 = -170 + 84*np.log(2) - 16*np.log(2)**2 + np.pi**2/3 * (10 - 4*np.log(2)) + 8*zeta(3)
    d1 = 88 - 40*np.log(2) + 8*np.log(2)**2 - 4/3*np.pi**2
    d2 = -20 + 8*np.log(2)
    d3 = 8/3
    d = [d0, d1, d2, d3]

    f1 = 2.910
    f2 = 78.35
    f3 = 1837
    f = [f1, f2, f3]

    if kappa < 25:

        sum_c_term = 0

        for i in range(len(c)):
            sum_c_term += c[i] * (kappa - 2)**(i + 1)

        return np.pi / 12 * (kappa - 2)**4 / (1 + sum_c_term)
    
    elif kappa >= 25:

        sum_d_term = 0
        sum_f_term = 0 

        for i in range(len(d)):
            sum_d_term += d[i] * np.log(kappa)**i 
        
        for i in range(len(f)):
            sum_f_term += f[i] * kappa**-(i + 1)

        return kappa * sum_d_term / (1 - sum_f_term)

# ----------------------------------------------------------------------------------------------------
def integrand_beta(kappa, Gmm, z):

    eps = kappa * const.m_e * const.c**2 / (2 * Gmm)

    return photon_density(eps, z) * phi(kappa) / kappa**2

# ----------------------------------------------------------------------------------------------------
def compute_pairproduction_beta(A, Z, Gmm, z): # 1 / s 

    b = const.alpha * r0**2 * const.c * Z**2 * const.m_e / (A * const.m_p) / Gmm
    b = b * quad(integrand_beta, 2, 1e4, args = (Gmm, z))[0] * const.m_e * const.c**2 
    return b

# ----------------------------------------------------------------------------------------------------