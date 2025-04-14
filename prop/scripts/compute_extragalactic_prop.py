from compute_pairproduction_beta import compute_pairproduction_beta
from compute_photopion_beta import compute_photopion_beta
from scipy.integrate import odeint
import numpy as np 

RESULTS_DIR = "../results"

E_arr = np.logspace(17, 22, num = 100)

km_to_Mpc = 3.2408e-20

c = 2.99792458e10 # cm / s
H0 = 67.4 * km_to_Mpc # s^-1
mp = 0.938e9 # eV 
OmgLbd = 0.685
OmgM = 0.315

# ----------------------------------------------------------------------------------------------------
def H(z): # Hubble parameter

    return H0 * np.sqrt((1 + z)**3 * OmgM + OmgLbd)

# ----------------------------------------------------------------------------------------------------
def abs_dt_dz(z):

    return 1 / (1 + z) / H(z)

# ----------------------------------------------------------------------------------------------------
def dE_dz(E, z):

    return E[0] / (1 + z) + E[0] * abs_dt_dz(z) * (compute_pairproduction_beta(1, 1, E[0] / mp, z) + compute_photopion_beta(E[0] / mp, z))

# ----------------------------------------------------------------------------------------------------
def write_Eg_vs_E(zg):

    Eg = np.zeros_like(E_arr)
    z_arr = np.linspace(0, zg, num = 50)

    for iEg, E in enumerate(E_arr):
        Eg[iEg] = odeint(dE_dz, E, z_arr)[-1]

    if zg == int(zg):
        zg_str = f"{int(zg)}"
    else:
        zg_str = str(zg).replace(".", "_")

    np.savetxt(f"{RESULTS_DIR}/Eg_vs_E_xchecks_zg{zg_str}.dat", np.column_stack((E_arr, Eg)), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
def write_Eg_vs_zg(E):

    z_arr = np.linspace(0, 1, num = 1000)
    np.savetxt(f"{RESULTS_DIR}/Eg_vs_zg_xchecks_E{int(np.log10(E))}.dat", np.column_stack((z_arr, odeint(dE_dz, E, z_arr))), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for zg in [0.05, 0.5, 1, 2, 3]:
        write_Eg_vs_E(zg)

    for E in [1e17, 1e18, 1e19, 1e20, 1e21]:
        write_Eg_vs_zg(E)

# ----------------------------------------------------------------------------------------------------