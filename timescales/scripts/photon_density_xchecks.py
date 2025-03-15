from matplotlib import lines
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'legend.fontsize': 'large',
'legend.title_fontsize': 'large',
'axes.labelsize': 'x-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

FIGURES_DIR = "../figures"
REFERENCES_DIR = "../references"
RESULTS_DIR = "../results"

eV_to_GeV = 1e-9
GeV_to_eV = 1e9
pc_to_cm = 3.0857e18

c = 3.e10           # cm/s
h = 4.135667696e-15 # eV.s
kB = 8.6173303e-5   # eV/K

D = 3.7e6 * pc_to_cm # cm  
R = 225 * pc_to_cm   # cm 

T_IR = 3.5e-3 / kB    # K
T_OPT = 332.5e-3 / kB # K
NORM_IR = 0.15268288372409347
NORM_OPT = 2.810860522825165e-09

# ----------------------------------------------------------------------------------------------------
def photon_density(eps, norm, T): # Number of photons per unit volume and energy

    return 8 * np.pi / (h * c)**3 * eps**2 * norm / (np.exp(eps/(kB * T)) - 1) # cm^-3 eV^-1

# ----------------------------------------------------------------------------------------------------
def compute_photon_density_ref(photon_spectrum): # cm^-3 eV^-1

    data = np.loadtxt(f"{REFERENCES_DIR}/SED_Condorelli_{photon_spectrum}.dat")
    photon_density = 9/4 * (D/R)**2 / c * data[:,1] / (data[:,0] * eV_to_GeV)**2 / GeV_to_eV

    return data[:,0], photon_density

# ----------------------------------------------------------------------------------------------------
def plot_photon_density_comparison():

    data_Condo_IR = compute_photon_density_ref('IR')
    data_Condo_OPT = compute_photon_density_ref('OPT')

    eps = np.logspace(-4, 1, num = 100)
    data_LAD_IR = photon_density(eps, NORM_IR, T_IR)
    data_LAD_OPT = photon_density(eps, NORM_OPT, T_OPT)

    Condo = lines.Line2D([], [], color = 'k', ls = '--', label = 'Condo+23')
    LAD = lines.Line2D([], [], color = 'k', ls = '-', label = 'LAD')
    lgnd = plt.legend(title = 'Results', handles = [Condo, LAD], frameon = True, loc = 'upper right')	
    plt.gca().add_artist(lgnd)

    plt.plot(data_Condo_IR[0], data_Condo_IR[1], ls = '--', color = '#FF5555')
    plt.plot(data_Condo_OPT[0], data_Condo_OPT[1], ls = '--', color = '#5555FF')
    plt.plot(eps, data_LAD_IR, color = '#FF5555', label = 'IR')
    plt.plot(eps, data_LAD_OPT, color = '#5555FF', label = 'OPT')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([1e-2, 1e8])
    plt.xlabel(r'Photon energy$\: \rm [eV]$')
    plt.ylabel(r'Photon density$\: \rm [cm^{-3} \, eV^{-1}]$')
    plt.legend(title = 'Radiation', loc = 'upper right', bbox_to_anchor = (1., 0.767))
    plt.savefig(f"{FIGURES_DIR}/photon_density_comparison.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/photon_density_comparison.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_photon_density_comparison()

# ----------------------------------------------------------------------------------------------------
