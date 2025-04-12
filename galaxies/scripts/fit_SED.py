from scipy.optimize import curve_fit
import numpy as np
import pandas as pd

DATA_DIR = "../data"
RESULTS_DIR = "../results"

erg_to_eV = 6.242e11 
eV_to_GeV = 1.e-9    
Jy_to_CGS = 1e-23 # erg·cm^-2·s^-1·Hz^-1
pc_to_cm = 3.0857e18 

D = 3.7e6 * pc_to_cm # cm  
R = 225 * pc_to_cm   # cm 

c = 3.e10           # cm/s
h = 4.135667696e-15 # eV·s
kB = 8.6173303e-5   # eV/K

# ----------------------------------------------------------------------------------------------------
def load_SED_measurements(galaxy):

    data = pd.read_csv(f"{DATA_DIR}/data_{galaxy}.csv")

    frequency = data["Frequency (Hz)"].to_numpy() 
    energy = frequency * h
    flux_density = data["Flux Density"].to_numpy() # Jy = 10^-26 W·m^-2·Hz^-1
    flux_density = flux_density * Jy_to_CGS * erg_to_eV * eV_to_GeV

    return energy, frequency*flux_density

# ----------------------------------------------------------------------------------------------------
def photon_density(E, Cdil, E0, T, sgm): # 1 / eV / cm^3

    return Cdil * 8 * np.pi * E**2 / (h * c)**3 / (np.exp(E / (kB * T)) - 1) * (E / E0)**sgm

# ----------------------------------------------------------------------------------------------------
def SED(E, Cdil, E0, T, sgm): 

    return 4/9 * (R/D)**2 * c * photon_density(E, Cdil, E0, T, sgm) * E**2 * eV_to_GeV

# ----------------------------------------------------------------------------------------------------
def fit_SED(galaxy):

    x_data, y_data = load_SED_measurements(galaxy)

    mask_energy = (x_data >= 1e-4) & (x_data <= 1e-1)
    mask_finite = np.isfinite(x_data) & np.isfinite(y_data)
    mask = mask_energy & mask_finite

    x_data = x_data[mask]
    y_data = y_data[mask]
    
    initial_guess = [0.1, 3.5e-3, 3.5e-3 / kB, 1.0]

    lower_bounds = [0, 0, 0, 0]
    upper_bounds = [np.inf, np.inf, np.inf, np.inf]

    return curve_fit(SED, x_data, y_data, p0 = initial_guess, bounds = (lower_bounds, upper_bounds))

# ----------------------------------------------------------------------------------------------------
def write_fitted_parameters(galaxy):

    header = 'Cdil' + '\t' + 'E0' + '\t' + 'T' + '\t' + 'sgm'
    np.savetxt(f"{RESULTS_DIR}/fitted_params_{galaxy}.dat", np.column_stack((fit_SED(galaxy)[0])), fmt = "%.15e", delimiter = "\t", header = header)

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    
    write_fitted_parameters('NGC253')

# ----------------------------------------------------------------------------------------------------