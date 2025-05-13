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

PARTICLES = ['1H', '4He', '14N', '28Si', '56Fe']
ZS = [1, 2, 7, 14, 26]

# ----------------------------------------------------------------------------------------------------
def iZ(Z):

    try:
        return ZS.index(Z)
    except ValueError:
        raise ValueError(f"Z ({Z}) not found in ZS.")

# ----------------------------------------------------------------------------------------------------
def injection_term(E, Z):

    return E**-2

# ----------------------------------------------------------------------------------------------------
def plot_cr_equilibrium_density(Z):

    data = np.loadtxt(f"{RESULTS_DIR}/cr_equilibrium_density_{PARTICLES[iZ(Z)]}.dat")

    plt.plot(np.log10(data[:,0]), data[:,0]**2 * data[:,1])
    plt.yscale('log')
    plt.xlabel(r'$\log_{10}(\rm Energy/eV)$')
    plt.ylabel(r'$E^{-2} n(E) \:$[arb. units]')
    plt.savefig(f"{FIGURES_DIR}/cr_equilibrium_density_{PARTICLES[iZ(Z)]}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/cr_equilibrium_density_{PARTICLES[iZ(Z)]}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
def plot_cr_escaping_emissivity(Z):

    data = np.loadtxt(f"{RESULTS_DIR}/cr_escaping_emissivity_{PARTICLES[iZ(Z)]}.dat")

    plt.plot(np.log10(data[:,0]), data[:,0]**2 * data[:,1], label = 'Escaping')
    plt.plot(np.log10(data[:,0]), data[:,0]**2 * injection_term(data[:,0], Z), label = 'Injected')
    plt.yscale('log')
    plt.xlabel(r'$\log_{10}(\rm Energy/eV)$')
    plt.ylabel(r'$E^{-2} Q_{\rm esc}(E) \:$[arb. units]')
    plt.legend()
    plt.savefig(f"{FIGURES_DIR}/cr_escaping_emissivity_{PARTICLES[iZ(Z)]}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/cr_escaping_emissivity_{PARTICLES[iZ(Z)]}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_cr_equilibrium_density(1)
    plot_cr_equilibrium_density(26)

    plot_cr_escaping_emissivity(1)
    plot_cr_escaping_emissivity(26)

# ----------------------------------------------------------------------------------------------------