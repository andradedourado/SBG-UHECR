from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd

plt.rcParams.update({'legend.fontsize': 'large',
'legend.title_fontsize': 'large',
'axes.labelsize': 'x-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

DATA_DIR = "../data"
FIGURES_DIR = "../figures"
RESULTS_DIR =  "../results"

gal_data = np.genfromtxt('starburst_galaxies.dat', dtype = None, encoding = None)

erg_to_eV = 6.242e11
eV_to_GeV = 1.e-9
Jy_to_CGS = 1e-23 # erg·cm^-2·s^-1·Hz^-1    
pc_to_cm = 3.0857e18 

R = 225 * pc_to_cm  # cm 

c = 3.e10           # cm/s
h = 4.135667696e-15 # eV·s
kB = 8.6173303e-5   # eV/K

# ----------------------------------------------------------------------------------------------------
def get_spectrum_color(photon_spectrum):
    
    if photon_spectrum == 'IR':
        return '#FF5555'
    
    elif photon_spectrum == 'OPT':
        return '#5555FF'

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
def plot_photon_spectrum_measurements(galaxy):

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

    plt.scatter(energy, frequency*flux_density, s = 25, marker = 'o', color = 'gray', edgecolor = 'black')

    return igal

# ----------------------------------------------------------------------------------------------------
def plot_fitted_SED(galaxy):

    E_IR = np.logspace(-4, -1, num = 100)
    E_OPT = np.logspace(-2, 1, num = 100)

    igal = plot_photon_spectrum_measurements(galaxy)

    for photon_spectrum in ('IR', 'OPT'):
        if photon_spectrum == 'IR':
            E = E_IR
        elif photon_spectrum == 'OPT':
            E = E_OPT
        fitted_params = np.loadtxt(f"{RESULTS_DIR}/fitted_params_{igal:02d}_{galaxy}_{photon_spectrum}.dat")
        plt.plot(E, SED(E, fitted_params[0], fitted_params[1], fitted_params[2], galaxy), color = get_spectrum_color(photon_spectrum), label = f'{photon_spectrum}')

    plt.gca().add_artist(AnchoredText(f'{galaxy}', loc = 'upper left', frameon = False, prop = {'fontsize': 'x-large'}))

    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1.e-6, 1.e2])
    plt.ylim([1.e-10, 2.e-4])
    plt.xlabel(r'Energy$\: \rm [eV]$')
    plt.ylabel(r'$E^2 F(E) \: \rm [GeV \, cm^{-2} \, s^{-1}]$')
    plt.legend(loc = 'lower left')
    plt.savefig(f"{FIGURES_DIR}/fitted_SED_{igal:02d}_{galaxy}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/fitted_SED_{igal:02d}_{galaxy}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for galaxy in gal_data:
        if galaxy[2] == 'NGC4038/9':
            plot_fitted_SED('NGC4038')
            plot_fitted_SED('NGC4039')
        else:
            plot_fitted_SED(galaxy[2])

# ----------------------------------------------------------------------------------------------------