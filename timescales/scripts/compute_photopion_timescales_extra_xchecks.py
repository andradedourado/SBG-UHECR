from scipy.integrate import simps
import numpy as np 

RESULTS_DIR = "../results"
XSECTIONS_DIR = "../../cross-sections/scripts"
PARTICLES = ['1H', '4He', '14N', '28Si', '56Fe']
ZS = [1, 2, 7, 14, 26]

mbarn_to_m2 = 1.e-31
GeV_to_eV = 1.e9
MeV_to_eV = 1.e6     
s_to_yr = 1 / (60 * 60 * 24 * 365.25) 

c =  299792458       # m/s
hbar = 6.5821220e-16 # eV.s
kB = 8.6173303e-5    # eV/K
mp = 0.9383e9        # eV 

T = 2.7 # K

# ----------------------------------------------------------------------------------------------------
def iZ(Z):

    try:
        return ZS.index(Z)
    except ValueError:
        raise ValueError(f"Z ({Z}) not found in ZS.")

# ----------------------------------------------------------------------------------------------------
def I(eps, Gmm):

    return -kB / (np.pi**2 * (hbar * c)**3) * T * np.log(1. - np.exp(-(eps)/(2 * Gmm * kB * T)))

# ----------------------------------------------------------------------------------------------------
def inelasticity(x): 

    Y_inf = 0.47
    xb = 6e9 # eV
    dlt = 0.33
    s = 0.15

    return Y_inf * (x / xb)**dlt / (1 + (x / xb)**(dlt / s))**s

# ----------------------------------------------------------------------------------------------------
def get_eps_and_cross_sections(A, Z):
    
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
    cross_section = A * np.array(cross_section) * mbarn_to_m2
    eps = (s - mp**2) / (2*mp)

    return eps, cross_section

# ----------------------------------------------------------------------------------------------------
def compute_timescales(A, Z, Gmms):

    timescales = []

    eps, cross_section = get_eps_and_cross_sections(A, Z)

    for Gmm in Gmms:
        integrand_interaction_rate = c / (2 * Gmm**2) * eps * cross_section * inelasticity(eps) * I(eps, Gmm)
        interaction_rate = simps(integrand_interaction_rate, eps)
        timescales.append(interaction_rate**-1 * s_to_yr)

    return np.array(timescales)

# ----------------------------------------------------------------------------------------------------
def write_timescales(A, Z):

    Gmms = np.logspace(8, 12, num = 100) / A
    E = Gmms * A * mp
    timescales = compute_timescales(A, Z, Gmms)   
    np.savetxt(f"{RESULTS_DIR}/timescales_photopion_{PARTICLES[iZ(Z)]}_CMB.dat", np.column_stack((E, timescales)), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_timescales(1, 1)

# ----------------------------------------------------------------------------------------------------