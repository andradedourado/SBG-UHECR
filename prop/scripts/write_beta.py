from compute_pairproduction_beta import compute_pairproduction_beta
from compute_photopion_beta import compute_photopion_beta
import numpy as np

RESULTS_DIR = "../results"

E_arr = np.logspace(17, 25, num = 100)

mp = 0.938e9 # eV

# ----------------------------------------------------------------------------------------------------
def write_beta_xchecks():

    beta = np.zeros_like(E_arr)
    
    for iE, E in enumerate(E_arr):
        beta[iE] = compute_pairproduction_beta(E / mp) + compute_photopion_beta(E / mp)

    np.savetxt(f"{RESULTS_DIR}/beta.dat", np.column_stack((E_arr, beta)), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_beta_xchecks()

# ----------------------------------------------------------------------------------------------------