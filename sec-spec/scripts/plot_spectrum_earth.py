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

pc_to_cm = 3.086e18
Mpc_to_cm = 3.086e24

D = 3.7 * Mpc_to_cm
R = 200 * pc_to_cm
V = 4 / 3 * np.pi * R**3

# ----------------------------------------------------------------------------------------------------
def compute_spectrum_earth(spec):

    return spec * V / (4 * np.pi * D**2)

# ----------------------------------------------------------------------------------------------------
def plot_spectrum_earth(l):
    
    if l == 'gmm':

        data_KAB06_gmm = np.loadtxt(f"{RESULTS_DIR}/KAB06_spectrum_gmm.dat")
        data_KA08_gmm = np.loadtxt(f"{RESULTS_DIR}/KA08_spectrum_gmm.dat")

        plt.plot(data_KAB06_gmm[:,0], compute_spectrum_earth(data_KAB06_gmm[:,1]) * data_KAB06_gmm[:,0]**2, label = r'$\gamma \: (pp)$')  
        plt.plot(data_KA08_gmm[:,0], compute_spectrum_earth(data_KA08_gmm[:,1]) * data_KA08_gmm[:,0]**2, label = r'$\gamma \: (p\gamma)$')

    if l == 'nu':

        data_KAB06_nu_mu = np.loadtxt(f"{RESULTS_DIR}/KAB06_spectrum_nu_mu.dat")
        data_KA08_anu_mu = np.loadtxt(f"{RESULTS_DIR}/KA08_spectrum_anu_mu.dat")
        data_KA08_nu_mu = np.loadtxt(f"{RESULTS_DIR}/KA08_spectrum_nu_mu.dat")
        data_KA08_nu_e = np.loadtxt(f"{RESULTS_DIR}/KA08_spectrum_nu_e.dat")
        data_KA08_anu_e = np.loadtxt(f"{RESULTS_DIR}/KA08_spectrum_anu_e.dat")

        plt.plot(data_KAB06_nu_mu[:,0], compute_spectrum_earth(data_KAB06_nu_mu[:,1]) * data_KAB06_nu_mu[:,0]**2, label = r'$\nu_\mu \: (pp)$')  
        plt.plot(data_KA08_anu_mu[:,0], compute_spectrum_earth(data_KA08_anu_mu[:,1]) * data_KA08_anu_mu[:,0]**2, label = r'$\bar{\nu}_\mu \: (p\gamma)$')
        plt.plot(data_KA08_nu_mu[:,0], compute_spectrum_earth(data_KA08_nu_mu[:,1]) * data_KA08_nu_mu[:,0]**2, label = r'$\nu_\mu \: (p\gamma)$')
        plt.plot(data_KA08_nu_e[:,0], compute_spectrum_earth(data_KA08_nu_e[:,1]) * data_KA08_nu_e[:,0]**2, label = r'$\nu_e \: (p\gamma)$')
        plt.plot(data_KA08_anu_e[:,0], compute_spectrum_earth(data_KA08_anu_e[:,1]) * data_KA08_anu_e[:,0]**2, label = r'$\bar{\nu}_e \: (p\gamma)$')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Energy$\: \rm [eV]$')
    plt.ylabel(r'Spectrum$\: \rm [eV \, cm^{-2} \, s^{-1}]$')
    plt.legend()
    plt.savefig(f"{FIGURES_DIR}/spectrum_{l}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/spectrum_{l}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_spectrum_earth('gmm') # No absorptions
    plot_spectrum_earth('nu')

# ----------------------------------------------------------------------------------------------------