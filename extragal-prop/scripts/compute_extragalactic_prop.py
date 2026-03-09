from compute_pairproduction_beta import compute_pairproduction_beta
from compute_photopion_beta import compute_photopion_beta
from scipy.integrate import odeint, simps
from scipy.interpolate import interp1d
import numpy as np 

GALAXIES_RESULTS_DIR = "../../galaxies/results"
LEAKY_BOX_RESULTS_DIR = "../../leaky-box/results"
RESULTS_DIR = "../results"

data = np.loadtxt(f"{GALAXIES_RESULTS_DIR}/lum_dist_to_z.dat")

E_arr = np.logspace(17, 22, num = 100)

Mpc_to_cm = 3.086e24
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

    return E / (1 + z) + abs_dt_dz(z) * E * (1 + z)**3 * (compute_pairproduction_beta(E / mp * (1 + z)) + compute_photopion_beta(E / mp * (1 + z))) 

# ----------------------------------------------------------------------------------------------------
def write_Eg_vs_E(zg, idata):

    Eg = np.zeros_like(E_arr)
    z_arr = np.linspace(0, zg, num = 50)

    for iEg, E in enumerate(E_arr):
        Eg[iEg] = odeint(dE_dz, E, z_arr)[-1]

    np.savetxt(f"{RESULTS_DIR}/Eg_vs_E_{idata+2:02d}.dat", np.column_stack((E_arr, Eg)), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
def write_dEg_dE(idata):

    data = np.loadtxt(f"{RESULTS_DIR}/Eg_vs_E_{idata+2:02d}.dat")
    dEg_dE = data[:,1] * np.gradient(np.log(data[:,1]), data[:,0])

    np.savetxt(f"{RESULTS_DIR}/dEg_dE_{idata+2:02d}.dat", np.column_stack((E_arr, dEg_dE)), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
def injection_term(Eg):

    esc_emis_data = np.loadtxt(f"{LEAKY_BOX_RESULTS_DIR}/cr_escaping_emissivity_1H.dat")
    Q_esc = interp1d(esc_emis_data[:,0], esc_emis_data[:,1], kind = 'linear', fill_value = 'extrapolate')
    return Q_esc(Eg)

# ----------------------------------------------------------------------------------------------------
def write_single_source_solution(zg, idata): 

    Eg = np.loadtxt(f"{RESULTS_DIR}/Eg_vs_E_{idata+2:02d}.dat")[:,1]
    dEg_dE = np.loadtxt(f"{RESULTS_DIR}/dEg_dE_{idata+2:02d}.dat")[:,1]

    Q0 = injection_term(Eg)

    n = (1 + zg) * Q0 / (4 * np.pi * c * (data[idata, 0] * Mpc_to_cm)**2) * dEg_dE

    np.savetxt(f"{RESULTS_DIR}/uhecr_density_{idata+2:02d}_injSBG.dat", np.column_stack((E_arr, n)), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for idata, zg in enumerate(data[:,1]):
        # write_Eg_vs_E(zg, idata)
        # write_dEg_dE(idata)
        write_single_source_solution(zg, idata)

# ----------------------------------------------------------------------------------------------------