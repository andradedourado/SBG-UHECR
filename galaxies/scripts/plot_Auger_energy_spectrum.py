import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'legend.fontsize': 'large',
'legend.title_fontsize': 'large',
'axes.labelsize': 'x-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

DATA_DIR = "../data"
FIGURES_DIR = "../figures"
RESULTS_DIR = "../results"

# ----------------------------------------------------------------------------------------------------
def plot_Auger_energy_spectrum():

    data = np.loadtxt(f"{DATA_DIR}/CombinedSpectrum_ICRC2019.txt")

    threshold = 38e18 # eV
    mask = 10**data[:,0] > threshold

    E = 10**data[:,0]
    J = data[:,1]
    J_lower = data[:,2]
    J_upper = data[:,3]

    J_err = [J - J_lower, J_upper - J] 
    J_err_scaled = [E**3 * (J_err[0]), E**3 * (J_err[1])]

    plt.axvspan(2e16, threshold, color = 'gray', alpha = 0.2)
    plt.errorbar(E, E**3 * J, yerr = J_err_scaled, color = 'k', linestyle = 'none', marker = '.', label = 'Auger (ICRC 2019)')
    plt.errorbar(E[mask], E[mask]**3 * J[mask] * 0.1, yerr = J_err_scaled[0][mask] * 0.1, color = 'gray', linestyle = 'none', marker = '.')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([2e16, 2e20])
    plt.xlabel(r'Energy$\: \rm [eV]$')
    plt.ylabel(r'$E^3 \times {\rm Intensity} \: \rm [eV^2 \, km^{-2} \, sr^{-1} \, yr^{-1}]$')
    plt.legend(frameon = False)
    plt.savefig(f"{FIGURES_DIR}/Auger_energy_spectrum.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/Auger_energy_spectrum.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_Auger_energy_spectrum()

# ----------------------------------------------------------------------------------------------------