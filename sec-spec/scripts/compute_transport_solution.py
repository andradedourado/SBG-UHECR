from scipy.integrate import quad
from scipy.interpolate import interp1d
import numpy as np

RESULTS_DIR = "../results"
TIMESCALES_DIR = "../../timescales/results"

erg_to_eV = 6.242e11
pc_to_cm = 3.086e18
yr_to_s = 60 * 60 * 24 * 365.25

R = 200 * pc_to_cm
R_SN = 0.05 * yr_to_s 
V = 4 / 3 * np.pi * R**3 # Volume 

Gmm = 2
# Gmm = 1.22
# Rcut = 10**18.72 # V

# ----------------------------------------------------------------------------------------------------
def injection_term_integrand(E, Z):

    return E**(-Gmm + 1) # * np.exp(-E / (Z * Rcut))

# ----------------------------------------------------------------------------------------------------
def injection_term(E, Z): # The rate of injection of particles per unit volume per unit time, Q(E)

    xi_CR = 0.1 
    E_SN = 1e51 * erg_to_eV 

    A = xi_CR * E_SN / quad(injection_term_integrand, 1e9, 1e21, args = (Z))[0] 
    return A * (R_SN / V) * E**-Gmm # * np.exp(-E / (Z * Rcut)) 

# ----------------------------------------------------------------------------------------------------
def get_total_timescale(E): # \tau
    
    files = ["timescales_advection.dat", "timescales_diff_1H.dat", "timescales_photopion_1H.dat", "timescales_pairproduction_1H.dat"]
    
    energy_conversion = {"timescales_advection.dat": 1.e18, "timescales_diff_1H.dat": 1.e18, "timescales_photopion_1H.dat": 1.0, "timescales_pairproduction_1H.dat": 1.0}

    interpolated_tau = []

    for file in files:
        data = np.loadtxt(f"{TIMESCALES_DIR}/{file}")
        E_table = data[:,0] * energy_conversion[file]
        tau_table = data[:,1] * yr_to_s  

        f_interp = interp1d(E_table, tau_table, kind = 'cubic', bounds_error = False, fill_value = np.inf)
        interpolated_tau.append(f_interp(E))

    total_tau_inv = sum(1 / tau for tau in interpolated_tau)
    return np.where(total_tau_inv > 0, 1 / total_tau_inv, np.inf)  

# ----------------------------------------------------------------------------------------------------
def write_transport_equation_solution(): # Number of particles per unit volume, n(E)

    E = np.logspace(17, 20, num = 500)
    np.savetxt(f"{RESULTS_DIR}/transport_sol_1H.dat", np.column_stack((E, injection_term(E, 1) / get_total_timescale(E))), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_transport_equation_solution()

# ----------------------------------------------------------------------------------------------------
