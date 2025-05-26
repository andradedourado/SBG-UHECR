from matplotlib import pyplot as plt
from scipy.integrate import quad
import numpy as np 

plt.rcParams.update({'legend.fontsize': 'large',
'legend.title_fontsize': 'large',
'axes.labelsize': 'x-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

FIGURES_DIR = "../figures"

PARTICLES = ['1H', '4He', '14N', '28Si', '56Fe']
PARTICLES_LEGEND = [r'$^1\mathrm{H}$', r'$^4\mathrm{He}$', r'$^{14}\mathrm{N}$', r'$^{28}\mathrm{Si}$', r'$^{56}\mathrm{Fe}$']
ZS = [1, 2, 7, 14, 26]

eV_to_erg = 1.60218e-12

Emin = 10**17.8 # eV
L0 = 5e44 # erg Mpc^-3 yr^-1
Gmm = -1.47
Rcut = 10**18.19 # V

E0 = 1e18

# ----------------------------------------------------------------------------------------------------
def iZ(Z):

    try:
        return ZS.index(Z)
    except ValueError:
        raise ValueError(f"Z ({Z}) not found in ZS.")

# ----------------------------------------------------------------------------------------------------
def integrand_generation_rate(E, Z):

    if E <= Z * Rcut: 
        return E * (E / E0)**-Gmm
    elif E > Z * Rcut: 
        return E * (E / E0)**-Gmm * np.exp(1 - E / (Z * Rcut)) 
    
# ----------------------------------------------------------------------------------------------------
def compute_generation_rate(E, Z):

    mask_low = E <= Z * Rcut
    mask_high = ~mask_low

    Q_A = np.zeros_like(E)

    Q_A[mask_low] = (E[mask_low] / E0)**-Gmm
    Q_A[mask_high] = (E[mask_high] / E0)**-Gmm * np.exp(1 - E[mask_high] / (Z * Rcut))

    return Q_A * [0.0, 0.245, 0.681, 0.049, 0.025][iZ(Z)] * L0 / (quad(integrand_generation_rate, Emin, 1e23, args = (Z))[0] * eV_to_erg**2)

# ----------------------------------------------------------------------------------------------------
def plot_generation_rate():

    E = np.logspace(17.8, 20.2, num = 100)

    for Z in ZS[1:]:
        plt.plot(np.log10(E), compute_generation_rate(E, Z), label = f'{PARTICLES_LEGEND[iZ(Z)]}')

    plt.yscale('log')
    plt.xlim([17.8, 20.2])
    plt.ylim(bottom = 1e26)
    plt.xlabel(r'$\log_{10}(\rm Energy/eV)$')
    plt.ylabel(r'$\tilde{Q}(E) \: \rm [erg^{-1} \, Mpc^{-3} \, yr^{-1}]$')
    plt.legend()
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_generation_rate()

# ----------------------------------------------------------------------------------------------------