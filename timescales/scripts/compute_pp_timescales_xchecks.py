from scipy.integrate import simps
from scipy.special import zeta
import matplotlib.pyplot as plt
import numpy as np 
import scipy.constants as const

REFERENCES_DIR = "../references"
RESULTS_DIR = "../results"
XSECTIONS_DIR = "../../cross-sections/scripts"

cm2_to_m2 = 1e-4 
eV_to_GeV = 1e-9
eV_to_J = 1.60218e-19
GeV_to_eV = 1.e9
pc_to_cm = 3.0857e18
s_to_yr = 1 / (60 * 60 * 24 * 365.25)

c =  299792458 # m/s
m_p = 0.9383e9  # eV
r0 = (const.e**2) / (4 * const.pi * const.epsilon_0 * const.m_e * const.c**2)

D = 3.7e6 * pc_to_cm # cm  
R = 225 * pc_to_cm   # cm 

# ----------------------------------------------------------------------------------------------------
def phi(kappa):

    phi = np.zeros_like(kappa)

    mask_1 = kappa < 25
    mask_2 = kappa >= 25

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

    sum_c_term = np.zeros_like(kappa)
    sum_d_term = np.zeros_like(kappa)
    sum_f_term = np.zeros_like(kappa)

    for i in range(len(c)):
        sum_c_term[mask_1] += c[i] * (kappa[mask_1] - 2)**(i + 1)
    
    for i in range(len(d)):
        sum_d_term[mask_2] += d[i] * np.log(kappa[mask_2])**i 
        
    for i in range(len(f)):
        sum_f_term[mask_2] += f[i] * kappa[mask_2]**-(i + 1)

    phi[mask_1] = np.pi / 12 * (kappa[mask_1] - 2)**4 / (1 + sum_c_term[mask_1])
    phi[mask_2] = kappa[mask_2] * sum_d_term[mask_2] / (1 - sum_f_term[mask_2])

    return phi

# ----------------------------------------------------------------------------------------------------
def compute_pairproduction_timescales(A, Z, Gmms):

    data_IR = np.loadtxt(f"{REFERENCES_DIR}/SED_Condorelli_IR.dat")
    data_OPT = np.loadtxt(f"{REFERENCES_DIR}/SED_Condorelli_OPT.dat")

    eps_IR = data_IR[:,0] * eV_to_J
    eps_OPT = data_OPT[:,0] * eV_to_J
    
    photon_density_IR = 9/4 * (D/R)**2 / c * data_IR[:,1] / (data_IR[:,0] * eV_to_GeV)**2 / GeV_to_eV / eV_to_J / cm2_to_m2
    photon_density_OPT = 9/4 * (D/R)**2 / c * data_OPT[:,1] / (data_OPT[:,0] * eV_to_GeV)**2 / GeV_to_eV / eV_to_J / cm2_to_m2
    
    timescales = []

    for Gmm in Gmms:

        kappa_IR = 2 * Gmm * eps_IR / (const.m_e * const.c**2)
        kappa_OPT = 2 * Gmm * eps_OPT / (const.m_e * const.c**2)

        integrand_IR = photon_density_IR * phi(kappa_IR) / kappa_IR**2
        integrand_OPT = photon_density_OPT * phi(kappa_OPT) / kappa_OPT**2

        mask_IR = (kappa_IR >= 2) & (kappa_IR < 1e4)
        mask_OPT = (kappa_OPT >= 2) & (kappa_OPT < 1e4)

        interaction_rate = const.alpha * r0**2 * const.c * Z**2 * const.m_e / (A * const.m_p) / Gmm        

        if np.count_nonzero(photon_density_OPT[mask_OPT] != 0) >= 2:
            interaction_rate = interaction_rate * (simps(integrand_IR[mask_IR], kappa_IR[mask_IR]) + simps(integrand_OPT[mask_OPT], kappa_OPT[mask_OPT])) * const.m_e * const.c**2
            timescales.append(interaction_rate**-1 * s_to_yr)
        else:
            interaction_rate = interaction_rate * simps(integrand_IR[mask_IR]) * const.m_e * const.c**2
            timescales.append(interaction_rate**-1 * s_to_yr)

    return np.array(timescales)

# ----------------------------------------------------------------------------------------------------
def write_pairproduction_timescales(A, Z):

    Gmms = np.logspace(8, 12, num = 100) / A
    E = Gmms * A * m_p 
    timescales = compute_pairproduction_timescales(A, Z, Gmms)   
    np.savetxt(f"{RESULTS_DIR}/timescales_pp_xchecks_1H.dat", np.column_stack((E, timescales)), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_pairproduction_timescales(1, 1)

# ----------------------------------------------------------------------------------------------------