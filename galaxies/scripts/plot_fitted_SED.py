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

erg_to_eV = 6.242e11
eV_to_GeV = 1.e-9
Jy_to_CGS = 1e-23 # erg·cm^-2·s^-1·Hz^-1    
pc_to_cm = 3.0857e18 

D = 3.7e6 * pc_to_cm # cm  
R = 225 * pc_to_cm  # cm 

c = 3.e10           # cm/s
h = 4.135667696e-15 # eV·s
kB = 8.6173303e-5   # eV/K

# ----------------------------------------------------------------------------------------------------
def photon_density(E, Cdil, E0, T, sgm): # 1 / eV / cm^3

    return Cdil * 8 * np.pi * E**2 / (h * c)**3 / (np.exp(E / (kB * T)) - 1) * (E / E0)**sgm

# ----------------------------------------------------------------------------------------------------
def SED(E, Cdil, E0, T, sgm): 

    return 4/9 * (R/D)**2 * c * photon_density(E, Cdil, E0, T, sgm) * E**2 * eV_to_GeV

# ----------------------------------------------------------------------------------------------------
def plot_photon_spectrum_measurements(galaxy):

    data = pd.read_csv(f"{DATA_DIR}/data_{galaxy}.csv")

    frequency = data["Frequency (Hz)"].to_numpy() 
    energy = frequency * h
    flux_density = data["Flux Density"].to_numpy() # Jy = 10^-26 W·m^-2·Hz^-1
    flux_density = flux_density * Jy_to_CGS * erg_to_eV * eV_to_GeV

    plt.scatter(energy, frequency*flux_density, s = 25, marker = 'o', color = 'gray', edgecolor = 'black')

# ----------------------------------------------------------------------------------------------------
def plot_fitted_SED_IR(galaxy):

    E = np.logspace(-4, -1, num = 100)

    tab10 = plt.get_cmap('tab10')
    tab10_colors = [tab10(i) for i in range(10)]

    plot_photon_spectrum_measurements(galaxy)

    fitted_params = np.loadtxt(f"{RESULTS_DIR}/fitted_params_{galaxy}.dat")
    plt.plot(E, SED(E, fitted_params[0], fitted_params[1], fitted_params[2], fitted_params[3]), color = tab10_colors[3], label = 'LAD (GB)')
    
    data = np.loadtxt(f"../../timescales/results/photon_spectrum_IR.dat")
    E, flux = data[:,0], data[:,1]  
    plt.plot(E * erg_to_eV, E**2 * flux * erg_to_eV * eV_to_GeV, color = tab10_colors[1], label = 'LAD (BB)')

    data = np.loadtxt(f"../../timescales/references/SED_Condorelli_IR.dat")
    plt.plot(data[:,0], data[:,1], color = tab10_colors[8], label = 'Condo+23')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1.e-6, 1.e2])
    plt.ylim([1.e-10, 2.e-4])
    plt.xlabel(r'Energy$\: \rm [eV]$')
    plt.ylabel(r'$E^2 F(E) \: \rm [GeV \, cm^{-2} \, s^{-1}]$')
    plt.legend()
    plt.savefig(f"{FIGURES_DIR}/fitted_SED_IR.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/fitted_SED_IR.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    
    plot_fitted_SED_IR('NGC253')

# ----------------------------------------------------------------------------------------------------