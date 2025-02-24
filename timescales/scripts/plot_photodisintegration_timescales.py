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

PARTICLES_LEGEND = [r'$^{14}\mathrm{N}$', r'$^{28}\mathrm{Si}$', r'$^{56}\mathrm{Fe}$']
PARTICLES = ['14N', '28Si', '56Fe']
ZS = [7, 14, 26]

# ----------------------------------------------------------------------------------------------------
def iZ(Z):

    try:
        return ZS.index(Z)
    except ValueError:
        raise ValueError(f"Z ({Z}) not found in ZS.")
    
# ----------------------------------------------------------------------------------------------------
def plot_timescales():

    for Z in ZS:
        data = np.loadtxt(f"{RESULTS_DIR}/timescales_photodisintegration_{PARTICLES[iZ(Z)]}.dat")
        plt.plot(np.log10(data[:,0]), data[:,1], label = f'{PARTICLES_LEGEND[iZ(Z)]}')

    plt.yscale('log')
    plt.xlabel(r'$\log_{10}{({\rm Energy}/{\rm eV})}$')
    plt.ylabel(r'Timescales$\: \rm [yr]$')
    plt.legend(title = 'Nucleus')
    plt.savefig(f"{FIGURES_DIR}/photodisintegration_timescales.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/photodisintegration_timescales.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_timescales()

# ----------------------------------------------------------------------------------------------------