from scipy.integrate import simps
import numpy as np

XSECTIONS_DIR = "../../cross-sections/scripts"

mbarn_to_m2 = 1.e-31
GeV_to_eV = 1.e9

c =  299792458       # m/s
hbar = 6.5821220e-16 # eV.s
kB = 8.6173303e-5    # eV/K
mp = 0.9383e9        # eV 

# ----------------------------------------------------------------------------------------------------
def I(eps, Gmm, z, T = 2.7):

    return -kB / (np.pi**2 * (hbar * c)**3) * T * (1 + z) * np.log(1. - np.exp(-(eps)/(2 * Gmm * kB * T * (1 + z)))) * (1 + z)**3

# ----------------------------------------------------------------------------------------------------
def inelasticity(x):

    Y_inf = 0.47
    xb = 6e9 # eV
    dlt = 0.33
    s = 0.15

    return Y_inf * (x / xb)**dlt / (1 + (x / xb)**(dlt / s))**s

# ----------------------------------------------------------------------------------------------------
def get_eps_and_cross_sections():

    s = []
    cross_section = []

    with open(f"{XSECTIONS_DIR}/xsecs_photopion_proton_sophia.txt", 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
                
            values = line.strip().split(',')
            s.append(float(values[0]))  
            cross_section.append(float(values[1]))

    s = np.array(s) * GeV_to_eV**2
    cross_section = np.array(cross_section) * mbarn_to_m2
    eps = (s - mp**2) / (2*mp)

    return eps, cross_section

# ----------------------------------------------------------------------------------------------------
def compute_photopion_beta(Gmm, z): # 1 / s

    eps, cross_section = get_eps_and_cross_sections()
    integrand_beta = c / (2 * Gmm**2) * eps * cross_section * inelasticity(eps) * I(eps, Gmm, z)
    return simps(integrand_beta, eps)

# ----------------------------------------------------------------------------------------------------