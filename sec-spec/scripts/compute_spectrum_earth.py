from scipy.interpolate import interp1d
import numpy as np 

GALAXIES_DIR = "../../galaxies/scripts"
RESULTS_DIR = "../results"
TIMESCALES_DIR = "../../timescales/results"

Mpc_to_cm = 3.086e24
yr_to_s = 60 * 60 * 24 * 365.25

# ----------------------------------------------------------------------------------------------------
def interpolate_spectrum(E, file):

    data = np.loadtxt(f"{RESULTS_DIR}/{file}")

    f_interp = interp1d(data[:,0], data[:,1], kind = 'cubic', bounds_error = False, fill_value = 0)

    return f_interp(E)

# ----------------------------------------------------------------------------------------------------
def compute_spectrum_earth(part):

    data_galaxies = np.genfromtxt(f"{GALAXIES_DIR}/starburst_galaxies.dat", dtype = None, encoding = None)

    E = np.logspace(13, 20, num = 100)
    spec = np.zeros_like(E)

    if part == 'gmm':
        files = ["KAB06_spectrum_gmm.dat", "KA08_spectrum_gmm.dat"]
    elif part == 'nu':
        files = ["KAB06_spectrum_nu_mu.dat", "KAB06_spectrum_nu_e.dat", "KA08_spectrum_anu_mu.dat", "KA08_spectrum_nu_mu.dat", "KA08_spectrum_nu_e.dat", "KA08_spectrum_anu_e.dat"]

    for igal in range(len(data_galaxies)):
        D = data_galaxies[igal][8] * Mpc_to_cm
        for file in files:
            spec += interpolate_spectrum(E, file) / (4 * np.pi * D**2)

    np.savetxt(f"{RESULTS_DIR}/spectrum_all_galaxies_{part}.dat", np.column_stack((E, spec)), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
def compute_proton_intensity(): # No extragalactic propagation
    
    data_galaxies = np.genfromtxt(f"{GALAXIES_DIR}/starburst_galaxies.dat", dtype = None, encoding = None)

    E = np.logspace(16, 20, num = 100)

    files = ["timescales_advection.dat", "timescales_diff_1H.dat"]

    interpolated_tau = []
    
    for file in files:
        data = np.loadtxt(f"{TIMESCALES_DIR}/{file}")
        E_table = data[:,0] * 1e18
        tau_table = data[:,1] * yr_to_s  
        f_interp = interp1d(E_table, tau_table, kind = 'cubic', bounds_error = False, fill_value = np.inf)
        interpolated_tau.append(f_interp(E))
    tau_inv_esc = sum(1 / tau for tau in interpolated_tau)
    tau_esc = 1 / tau_inv_esc

    Q_esc = np.loadtxt(f"{RESULTS_DIR}/transport_sol_1H.dat")[:,1] / tau_esc
    spec = np.zeros_like(Q_esc)

    for igal in range(len(data_galaxies)):
        D = data_galaxies[igal][8] * Mpc_to_cm
        spec += Q_esc / (4 * np.pi * D**2)

    np.savetxt(f"{RESULTS_DIR}/proton_flux.dat", np.column_stack((E, spec / (4 * np.pi))), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # compute_spectrum_earth('gmm')
    # compute_spectrum_earth('nu')
    
    compute_proton_intensity()

# ----------------------------------------------------------------------------------------------------