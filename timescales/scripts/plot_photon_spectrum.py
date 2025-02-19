from matplotlib import lines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({'legend.fontsize': 'large',
'legend.title_fontsize': 'large',
'axes.labelsize': 'x-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

FIGURES_DIR = "../figures"
RESULTS_DIR = "../results"
GALAXIES = ['NGC253', 'M82']

erg_to_eV = 6.242e11 # eV
eV_to_GeV = 1.e-9    # GeV
Jy_to_CGS = 1e-23    # erg·cm^-2·s^-1·Hz^-1

h = 4.135667696e-15 # eV·s

# ----------------------------------------------------------------------------------------------------
def get_marker_color_and_label(galaxy):

    if galaxy == 'NGC253':
        return 'o', 'gray', 'NGC 253'
    elif galaxy == 'M82':
        return 's', 'whitesmoke', 'M82'

# ----------------------------------------------------------------------------------------------------
def plot_photon_spectrum_measurements(galaxy):

    data = pd.read_csv(f"data_{galaxy}.csv")

    frequency = data["Frequency (Hz)"].to_numpy() 
    energy = frequency * h
    flux_density = data["Flux Density"].to_numpy() # Jy = 10^-26 W·m^-2·Hz^-1
    flux_density = flux_density * Jy_to_CGS * erg_to_eV * eV_to_GeV

    plt.scatter(energy, frequency*flux_density, s = 25, marker = get_marker_color_and_label(galaxy)[0], color = get_marker_color_and_label(galaxy)[1], edgecolor = 'black', label = get_marker_color_and_label(galaxy)[2])

# ----------------------------------------------------------------------------------------------------
def plot_photon_spectrum():

    data_IR = np.loadtxt(f"{RESULTS_DIR}/photon_spectrum_IR.dat") 
    data_OPT = np.loadtxt(f"{RESULTS_DIR}/photon_spectrum_OPT.dat")

    E_IR, intensity_IR = data_IR[:,0], data_IR[:,1]   
    E_OPT, intensity_OPT = data_OPT[:,0], data_OPT[:,1] 

    plt.figure()
    plt.plot(E_IR * erg_to_eV, E_IR**2 * intensity_IR * erg_to_eV * eV_to_GeV, color = '#FF5555')
    plt.plot(E_OPT * erg_to_eV, E_OPT**2 * intensity_OPT * erg_to_eV * eV_to_GeV, color = '#5555FF')
    plot_photon_spectrum_measurements('NGC253')
    plot_photon_spectrum_measurements('M82')

    IR = lines.Line2D([], [], color = '#FF5555', ls = '-', label = 'Far infrared')
    OPT = lines.Line2D([], [], color = '#5555FF', ls = '-', label = 'Optical')
    lgnd = plt.legend(title = 'Models', handles = [IR, OPT], frameon = True, loc = 'upper left')	
    plt.gca().add_artist(lgnd)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1.e-6, 1.e2])
    plt.ylim([1.e-10, 2.e-4])
    plt.xlabel(r'Energy$\: \rm [eV]$')
    plt.ylabel(r'$E^2 F(E) \: \rm [GeV \, cm^{-2} \, s^{-1}]$')
    plt.legend(title = 'Measurements', loc = 'upper left', bbox_to_anchor = (0., 0.769))
    plt.savefig(f"{FIGURES_DIR}/photon_spectrum.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/photon_spectrum.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_photon_spectrum()

# ----------------------------------------------------------------------------------------------------
