from scipy.interpolate import interp1d
import numpy as np

RESULTS_DIR = "../results"
SEC_SPEC_DIR = "../../sec-spec/results"
TIMESCALES_DIR = "../../timescales/results"

Mpc_to_cm = 3.086e24
yr_to_s = 60 * 60 * 24 * 365.25

# ----------------------------------------------------------------------------------------------------
def compute_escaping_particles_emissivity():
        
    E = np.logspace(16, 20, num = 500)

    files = ["timescales_advection.dat", "timescales_diff_1H.dat"]

    interpolated_tau = []
    
    for file in files:
        data = np.loadtxt(f"{TIMESCALES_DIR}/{file}")
        E_table = data[:,0] * 1e18
        tau_table = data[:,1] * yr_to_s  
        f_interp = interp1d(E_table, tau_table, kind = 'linear', bounds_error = False, fill_value = np.inf)
        interpolated_tau.append(f_interp(E))
    tau_inv_esc = sum(1 / tau for tau in interpolated_tau)
    tau_esc = 1 / tau_inv_esc

    return E, np.loadtxt(f"{SEC_SPEC_DIR}/transport_sol_1H.dat")[:,1] / tau_esc

# ----------------------------------------------------------------------------------------------------
def write_escaping_particles_emissivity():

    E, Q_esc = compute_escaping_particles_emissivity()

    np.savetxt(f"{RESULTS_DIR}/escaping_particles_emissivity.dat", np.column_stack((E, Q_esc)), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_escaping_particles_emissivity()

# ----------------------------------------------------------------------------------------------------
