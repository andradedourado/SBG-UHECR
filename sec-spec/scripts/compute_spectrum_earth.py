from scipy.interpolate import interp1d
import numpy as np 

GALAXIES_DIR = "../../galaxies/scripts"
RESULTS_DIR = "../results"

Mpc_to_cm = 3.086e24

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
if __name__ == '__main__':

    compute_spectrum_earth('gmm')
    # compute_spectrum_earth('nu')

# ----------------------------------------------------------------------------------------------------