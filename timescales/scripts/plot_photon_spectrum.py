import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'legend.fontsize': 'large',
'legend.title_fontsize': 'large',
'axes.labelsize': 'x-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

FIGURES_DIR = "../figures"
RESULTS_DIR = "../results"

erg_to_eV = 6.242e11 # eV
eV_to_GeV = 1.e-9    # GeV

# ----------------------------------------------------------------------------------------------------
def plot_photon_spectrum():

    data_IR = np.loadtxt(f"{RESULTS_DIR}/photon_spectrum_IR.dat") 
    data_OPT = np.loadtxt(f"{RESULTS_DIR}/photon_spectrum_OPT.dat")

    E_IR, intensity_IR = data_IR[:,0], data_IR[:,1]   
    E_OPT, intensity_OPT = data_OPT[:,0], data_OPT[:,1] 

    plt.figure()
    plt.plot(E_IR * erg_to_eV, E_IR**2 * intensity_IR * erg_to_eV * eV_to_GeV, label = 'Far infrared')
    plt.plot(E_OPT * erg_to_eV, E_OPT**2 * intensity_OPT * erg_to_eV * eV_to_GeV, label = 'Optical')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(bottom = 1.e-10)
    plt.xlabel(r'Energy$\: \rm [eV]$')
    plt.ylabel(r'$E^2 F(E) \: \rm [GeV \, cm^{-2} \, s^{-1}]$')
    plt.legend()
    plt.savefig(f"{FIGURES_DIR}/photon_spectrum.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/photon_spectrum.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_photon_spectrum()

# ----------------------------------------------------------------------------------------------------
