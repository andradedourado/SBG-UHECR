from compute_pairproduction_beta import compute_pairproduction_beta
from compute_photopion_beta import compute_photopion_beta
from scipy.integrate import odeint
# from scipy.interpolate import interp1d
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
def get_zg_arr(zg):

    if zg == int(zg):
        return f"{int(zg)}"
    else:
        return str(zg).replace(".", "_")

# ----------------------------------------------------------------------------------------------------
def H(z): # Hubble parameter

    return H0 * np.sqrt((1 + z)**3 * OmgM + OmgLbd)

# ----------------------------------------------------------------------------------------------------
def abs_dt_dz(z):

    return 1 / (1 + z) / H(z)

# ----------------------------------------------------------------------------------------------------
def dE_dz(E, z):

    # data = np.loadtxt(f"{RESULTS_DIR}/beta.dat")
    # interp_beta = interp1d(data[:,0], data[:,1], kind = 'linear', bounds_error = False, fill_value = np.nan)
    # return E / (1 + z) + abs_dt_dz(z) * E * (1 + z)**3 * interp_beta(E * (1 + z)) 

    return E / (1 + z) + abs_dt_dz(z) * E * (1 + z)**3 * (compute_pairproduction_beta(E / mp * (1 + z)) + compute_photopion_beta(E / mp * (1 + z))) 

# ----------------------------------------------------------------------------------------------------
def write_Eg_vs_E_xchecks(zg):

    Eg = np.zeros_like(E_arr)
    z_arr = np.linspace(0, zg, num = 50)

    for iEg, E in enumerate(E_arr):
        Eg[iEg] = odeint(dE_dz, E, z_arr)[-1]

    np.savetxt(f"{RESULTS_DIR}/Eg_vs_E_xchecks_zg{get_zg_arr(zg)}.dat", np.column_stack((E_arr, Eg)), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
def write_Eg_vs_z_xchecks(E):

    z_arr = np.linspace(0, 1, num = 100)
    np.savetxt(f"{RESULTS_DIR}/Eg_vs_z_xchecks_E{int(np.log10(E))}.dat", np.column_stack((z_arr, odeint(dE_dz, E, z_arr))), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for zg in [0.01, 0.05, 0.5, 1, 2, 3]:
        write_Eg_vs_E_xchecks(zg)

    for E in [1e17, 1e18, 1e19, 1e20, 1e21]:
        write_Eg_vs_z_xchecks(E)

# ----------------------------------------------------------------------------------------------------