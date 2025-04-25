from compute_pairproduction_beta import compute_pairproduction_beta
from compute_photopion_beta import compute_photopion_beta
from scipy.integrate import odeint, simps
from scipy.interpolate import interp1d
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
def derivative_b_pgamma_wrt_E_prime(z_arr): 

    b_pgamma = np.zeros((len(z_arr), len(E_arr))) 

    for iz, z in enumerate(z_arr): 
        for iE, E_prime in enumerate(E_arr):
            b_pgamma[iz, iE] = E_prime * (1 + z)**3 * (compute_pairproduction_beta(E_prime / mp) + compute_photopion_beta(E_prime / mp))

    return np.gradient(b_pgamma, E_arr, axis = 1)

# ----------------------------------------------------------------------------------------------------
def compute_dEg_dE(E, zg):

    z_arr = np.logspace(0, zg, num = 100)

    Eg = np.loadtxt(f"{RESULTS_DIR}/Eg_vs_E_xchecks_zg{get_zg_arr(zg)}.dat")[:,1] 
    
    interp_Eg = interp1d(E_arr, Eg, kind = 'linear', fill_value = 'extrapolate')
    Eg_val = interp_Eg(E)

    db_pgamma_grid = derivative_b_pgamma_wrt_E_prime(z_arr)
    db_pgamma_at_interp_Eg = np.array([interp1d(E_arr, db_pgamma_grid[iz], kind = 'linear', fill_value = 'extrapolate')((1 + z) * Eg_val) \
                            for iz, z in enumerate(z_arr)])
 
    integrand_dEg_dE = abs_dt_dz(z_arr) * db_pgamma_at_interp_Eg

    return (1 + zg) * np.exp(simps(integrand_dEg_dE, z_arr)) 

# ----------------------------------------------------------------------------------------------------
def write_dEg_dE(zg):

    dEg_dE = np.zeros_like(E_arr)

    for iE, E in enumerate(E_arr):
        dEg_dE[iE] = compute_dEg_dE(E, zg)

    np.savetxt(f"{RESULTS_DIR}/dEg_dE_zg{get_zg_arr(zg)}.dat", np.column_stack((E_arr, dEg_dE)), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
def injection_term(E):

    return E**-2

# ----------------------------------------------------------------------------------------------------
def write_single_source_solution(zg): 

    Eg = np.loadtxt(f"{RESULTS_DIR}/Eg_vs_E_xchecks_zg{get_zg_arr(zg)}.dat")[:,1]
    dEg_dE = np.loadtxt(f"{RESULTS_DIR}/dEg_dE_zg{get_zg_arr(zg)}.dat")[:,1]

    Q0 = injection_term(Eg)

    np.savetxt(f"{RESULTS_DIR}/single_source_solution_zg{get_zg_arr(zg)}.dat", np.column_stack((E_arr, Q0 * dEg_dE)), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # for zg in [0.01, 0.05, 0.5, 1, 2, 3]:
    #     write_Eg_vs_E_xchecks(zg)

    # for E in [1e17, 1e18, 1e19, 1e20, 1e21]:
    #     write_Eg_vs_z_xchecks(E)

    # for zg in [0.05, 0.5, 1, 2, 3]:
    #     write_dEg_dE(zg)

    for zg in [0.05, 0.5, 1, 2, 3]:
        write_single_source_solution(zg)    

# ----------------------------------------------------------------------------------------------------