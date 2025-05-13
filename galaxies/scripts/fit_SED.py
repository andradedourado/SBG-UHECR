from scipy.optimize import curve_fit
import numpy as np
import pandas as pd

DATA_DIR = "../data"
RESULTS_DIR = "../results"

gal_data = np.genfromtxt('starburst_galaxies.dat', dtype = None, encoding = None)

erg_to_eV = 6.242e11 
eV_to_GeV = 1.e-9    
Jy_to_CGS = 1e-23 # erg·cm^-2·s^-1·Hz^-1
pc_to_cm = 3.0857e18 

R = 225 * pc_to_cm   # cm 

c = 3.e10           # cm/s
h = 4.135667696e-15 # eV·s
kB = 8.6173303e-5   # eV/K

# ----------------------------------------------------------------------------------------------------
def get_initial_guess(photon_spectrum):

    if photon_spectrum == 'IR':
        return [0.1, 3.5e-3 / kB, 1.0]
    elif photon_spectrum == 'OPT':
        return [1e-09, 332.5e-3 / kB, 1.0]

# ----------------------------------------------------------------------------------------------------
def load_SED_measurements(galaxy):

    if galaxy in ('NGC4038', 'NGC4039'):
        igal = 36
        data = pd.read_csv(f"{DATA_DIR}/data_{igal}_{galaxy}.csv")
    
    else:
        for igal in range(len(gal_data)):
            if galaxy == gal_data[igal][2]:
                break
        igal = igal + 2
        data = pd.read_csv(f"{DATA_DIR}/data_{igal:02d}_{galaxy}.csv")

    frequency = data["Frequency (Hz)"].to_numpy() 
    energy = frequency * h
    flux_density = data["Flux Density"].to_numpy() # Jy = 10^-26 W·m^-2·Hz^-1
    flux_density = flux_density * Jy_to_CGS * erg_to_eV * eV_to_GeV

    return energy, frequency*flux_density, igal 

# ----------------------------------------------------------------------------------------------------
def photon_density(E, Cdil, T, sgm): # 1 / eV / cm^3

    return Cdil * 8 * np.pi * E**2 / (h * c)**3 / (np.exp(E / (kB * T)) - 1) * (E / (kB * T))**sgm

# ----------------------------------------------------------------------------------------------------
def SED(E, Cdil, T, sgm, galaxy): 

    if galaxy in ('NGC4038', 'NGC4039'):
        D = 23.7e6 * pc_to_cm
    else:
        D = float(gal_data[[g[2] for g in gal_data].index(galaxy)][8]) * 1e6 * pc_to_cm

    return 4/9 * (R/D)**2 * c * photon_density(E, Cdil, T, sgm) * E**2 * eV_to_GeV

# ----------------------------------------------------------------------------------------------------
def fit_SED(galaxy, photon_spectrum):

    x_data, y_data, igal = load_SED_measurements(galaxy)

    if photon_spectrum == 'IR':
        mask_energy = (x_data >= 1e-4) & (x_data <= 1e-1)
        initial_guess = get_initial_guess(photon_spectrum)
        fit_func = lambda E, Cdil, T, sgm: SED(E, Cdil, T, sgm, galaxy)
        bounds = ([0, 0, 0], [np.inf, np.inf, np.inf])
        p0 = initial_guess

    elif photon_spectrum == 'OPT':
        mask_energy = (x_data >= 1e-2) & (x_data <= 1e2)
        initial_guess = get_initial_guess(photon_spectrum)[:2]
        fit_func = lambda E, Cdil, T: SED(E, Cdil, T, 0, galaxy)
        bounds = ([0, 0], [np.inf, np.inf])
        p0 = initial_guess

    mask_finite = np.isfinite(x_data) & np.isfinite(y_data)
    mask = mask_energy & mask_finite

    x_data = x_data[mask]
    y_data = y_data[mask]

    if galaxy == 'NGC253':
        popt, _ = curve_fit(fit_func, x_data, y_data, p0 = p0, bounds = bounds)
    else:
        popt, _ = curve_fit(fit_func, x_data, y_data, p0 = p0, bounds = bounds, method = 'dogbox')

    if photon_spectrum == 'OPT':
        popt = np.append(popt, 0.0)

    return popt, igal

# ----------------------------------------------------------------------------------------------------
def write_fitted_parameters(galaxy, photon_spectrum):

    header = 'Cdil' + '\t' + 'T' + '\t' + 'sgm'

    if galaxy == 'NGC4038/9':
        np.savetxt(f"{RESULTS_DIR}/fitted_params_36_NGC4038_{photon_spectrum}.dat", fit_SED('NGC4038', photon_spectrum)[0].reshape(1, -1), fmt = "%.15e", delimiter = "\t", header = header)
        np.savetxt(f"{RESULTS_DIR}/fitted_params_36_NGC4039_{photon_spectrum}.dat", fit_SED('NGC4039', photon_spectrum)[0].reshape(1, -1), fmt = "%.15e", delimiter = "\t", header = header)
    else:
        igal = fit_SED(galaxy, photon_spectrum)[1]
        np.savetxt(f"{RESULTS_DIR}/fitted_params_{igal:02d}_{galaxy}_{photon_spectrum}.dat", fit_SED(galaxy, photon_spectrum)[0].reshape(1, -1), fmt = "%.15e", delimiter = "\t", header = header)

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for igal in range(len(gal_data)):
        write_fitted_parameters(gal_data[igal][2], 'IR')
        write_fitted_parameters(gal_data[igal][2], 'OPT')

# ----------------------------------------------------------------------------------------------------