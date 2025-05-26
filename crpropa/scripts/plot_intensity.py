import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'legend.fontsize': 'large',
'legend.title_fontsize': 'large',
'axes.labelsize': 'x-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

RESULTS_DIR = "../results"

PARTICLES = ['1H', '4He', '14N', '28Si', '56Fe']
PARTICLES_LEGEND = [r'$^1\mathrm{H}$', r'$^4\mathrm{He}$', r'$^{14}\mathrm{N}$', r'$^{28}\mathrm{Si}$', r'$^{56}\mathrm{Fe}$']
ZSS = [1, 2, 7, 14, 26]

# ----------------------------------------------------------------------------------------------------
def iZs(Zs):

    try:
        return ZSS.index(Zs)
    except ValueError:
        raise ValueError(f"Zs ({Zs}) not found in ZSS.")

# ----------------------------------------------------------------------------------------------------
def plot_intensity():

    for Zs in ZSS[1:]:        
        data = np.loadtxt(f"{RESULTS_DIR}/spec_{PARTICLES[iZs(Zs)]}.dat")
        plt.plot(np.log10(data[:,0]), data[:,0]**3 * data[:,1], label = f'{PARTICLES_LEGEND[iZs(Zs)]}')

    plt.yscale('log')
    plt.xlim([18, 21])
    plt.ylim([1e68, 1e71])
    plt.xlabel(r'$\log_{10}(\rm Energy/eV)$')
    plt.ylabel(r'$E^3 \times \rm Intensity \: [arb. units]$')
    plt.legend()
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_intensity()

# ----------------------------------------------------------------------------------------------------