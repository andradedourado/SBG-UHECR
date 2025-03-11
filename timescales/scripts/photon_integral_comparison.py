from scipy.integrate import quad
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'legend.fontsize': 'large',
'legend.title_fontsize': 'large',
'axes.labelsize': 'x-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

FIGURES_DIR = "../figures"

c =  299792458       # m/s
h = 4.135667696e-15  # eV.s 
hbar = 6.5821220e-16 # eV.s
kB = 8.6173303e-5    # eV/K

T_IR = 3.5e-3 / kB    # K
T_OPT = 332.5e-3 / kB # K
NORM_IR = 0.15268288372409347
NORM_OPT = 2.810860522825165e-09

# ----------------------------------------------------------------------------------------------------
def photon_density(eps): # Number of photons per unit volume per energy cubed for IR and OPT
    
    return 8 * np.pi / (h * c)**3 * (NORM_IR / (np.exp(eps / (kB * T_IR)) - 1) + NORM_OPT / (np.exp(eps / (kB * T_OPT)) - 1))

# ----------------------------------------------------------------------------------------------------
def I(eps_prime, Gmm):

    return -kB / (np.pi**2 * (hbar * c)**3) * (T_IR * NORM_IR * np.log(1. - np.exp(-(eps_prime)/(2 * Gmm * kB * T_IR))) + T_OPT * NORM_OPT * np.log(1. - np.exp(-(eps_prime)/(2 * Gmm * kB * T_OPT))))

# ----------------------------------------------------------------------------------------------------
def integrate_photon_density(eps_prime_arr, Gmm):

    return np.array([quad(photon_density, eps_prime / (2 * Gmm), np.inf)[0] for eps_prime in eps_prime_arr])

# ----------------------------------------------------------------------------------------------------
def plot_photon_integral_comparison():

    Gmm = 1.e10
    eps_prime_arr = np.logspace(0, 11, num = 100)

    plt.plot(eps_prime_arr, I(eps_prime_arr, Gmm), label = 'Analytical')
    plt.plot(eps_prime_arr, integrate_photon_density(eps_prime_arr, Gmm), label = 'Numerical')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r"$\epsilon' \: \rm [eV]$")
    plt.ylabel(r"Integral of photon density$\: \rm [m^{-3} eV^{-2}]$")
    plt.legend()
    plt.savefig(f"{FIGURES_DIR}/photon_integral_comparison.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/photon_integral_comparison.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_photon_integral_comparison()

# ----------------------------------------------------------------------------------------------------