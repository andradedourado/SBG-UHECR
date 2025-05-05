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
RESULTS_DIR = "../results"

gal_data = np.genfromtxt('starburst_galaxies.dat', dtype = None, encoding = None)

erg_to_eV = 6.242e11 
eV_to_GeV = 1.e-9    
Jy_to_CGS = 1e-23 # erg·cm^-2·s^-1·Hz^-1

h = 4.135667696e-15 # eV·s

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
def plot_SED(galaxy):

    E, SED, igal = load_SED_measurements(galaxy)

    plt.scatter(E, SED, s = 25, marker = 'o', color = 'gray', edgecolor = 'black', label = galaxy)

    plt.gca().add_artist(AnchoredText(f'{galaxy}', loc = 'upper left', frameon = False, prop = {'fontsize': 'x-large'}))

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Energy$\: \rm [eV]$')
    plt.ylabel(r'$E^2 F(E) \: \rm [GeV \, cm^{-2} \, s^{-1}]$') 
    plt.savefig(f"{FIGURES_DIR}/SED_{igal}_{galaxy}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/SED_{igal}_{galaxy}.png", bbox_inches = 'tight', dpi = 300)    
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for galaxy in gal_data:
        if galaxy[2] == 'NGC4038/9':
            plot_SED('NGC4038')
            plot_SED('NGC4039')
        else:
            plot_SED(galaxy[2])

# ----------------------------------------------------------------------------------------------------