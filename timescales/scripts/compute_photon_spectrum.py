from scipy.integrate import quad
import numpy as np

RESULTS_DIR = "../results"  

pc_to_cm = 3.0857e18 
eV_to_erg = 1.60218e-12 

D = 3.7e6 * pc_to_cm # cm  
R = 225 * pc_to_cm  # cm 

c = 3.e10          # cm/s
h = 6.62607015e-27 # erg·s
kB = 1.380649e-16  # erg/K

# ----------------------------------------------------------------------------------------------------
def photon_density(E, T): # Number of photons per unit volume

    return 8 * np.pi * E**3 / (h * c)**3 / (np.exp(E / (kB * T)) - 1)

# ----------------------------------------------------------------------------------------------------
def select_energy_range(photon_spectrum, nbins):

    if photon_spectrum == 'IR':
        return np.logspace(-5, -1, num = nbins) * eV_to_erg # erg
    elif photon_spectrum == 'OPT':
        return np.logspace(-3, 1, num = nbins) * eV_to_erg # erg
    else: 
        raise ValueError(f"Unknown photon spectrum: {photon_spectrum}. Valid options are 'IR' and 'OPT'.")

# ----------------------------------------------------------------------------------------------------
def select_blackbody_temperature(photon_spectrum): 

    if photon_spectrum == 'IR':
        return 3.5e-3 * eV_to_erg / kB # K  
    elif photon_spectrum == 'OPT':
        return 332.5e-3 * eV_to_erg / kB # K
    else: 
        raise ValueError(f"Unknown photon spectrum: {photon_spectrum}. Valid options are 'IR' and 'OPT'.")
    
# ----------------------------------------------------------------------------------------------------
def select_energy_density(photon_spectrum):

    if photon_spectrum == 'IR':
        return 1958 * eV_to_erg # erg cm^-3
    elif photon_spectrum == 'OPT':
        return 2936 * eV_to_erg # erg cm^-3
    else: 
        raise ValueError(f"Unknown photon spectrum: {photon_spectrum}. Valid options are 'IR' and 'OPT'.")

# ----------------------------------------------------------------------------------------------------
def normalize_photon_spectrum(photon_spectrum):
    
    if photon_spectrum == 'IR':
        integral = quad(photon_density, 10**-5 * eV_to_erg, 10**-1 * eV_to_erg, args = (select_blackbody_temperature(photon_spectrum)))[0]
    if photon_spectrum == 'OPT':
        integral = quad(photon_density, 10**-3 * eV_to_erg, 10**1 * eV_to_erg, args = (select_blackbody_temperature(photon_spectrum)))[0]

    return select_energy_density(photon_spectrum) / integral

# ----------------------------------------------------------------------------------------------------
def compute_photon_spectrum(photon_spectrum, nbins):

    E = select_energy_range(photon_spectrum, nbins) 
    T = select_blackbody_temperature(photon_spectrum)

    return E, 4/9 * (R/D)**2 * c * photon_density(E, T) * normalize_photon_spectrum(photon_spectrum) / E

# ----------------------------------------------------------------------------------------------------
def write_photon_spectrum(photon_spectrum, nbins = 100):

    E_values, spectrum_values = compute_photon_spectrum(photon_spectrum, nbins)

    with open(f"{RESULTS_DIR}/photon_spectrum_{photon_spectrum}.dat", 'w') as f:
        for E, spectrum in zip(E_values, spectrum_values):
            f.write(f'{E:.15e}\t{spectrum:.15e}\n')

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    
    write_photon_spectrum('IR')
    write_photon_spectrum('OPT')

# ----------------------------------------------------------------------------------------------------