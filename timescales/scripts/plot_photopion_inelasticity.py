import matplotlib.pyplot as plt
import numpy as np 

plt.rcParams.update({'legend.fontsize': 'large',
'legend.title_fontsize': 'large',
'axes.labelsize': 'x-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

FIGURES_DIR = "../figures"

# ----------------------------------------------------------------------------------------------------
def inelasticity(x): 

    Y_inf = 0.47
    xb = 6e9 # eV
    dlt = 0.33
    s = 0.15

    return Y_inf * (x / xb)**dlt / (1 + (x / xb)**(dlt / s))**s

# ----------------------------------------------------------------------------------------------------
def plot_photopion_inelasticity():

    x = np.logspace(8, 14, num = 100) # x = eps_prime

    plt.plot(x, inelasticity(x))
    plt.xscale('log')
    plt.xlim([1e8, 4e13])
    plt.ylim([0, 0.7])
    plt.xlabel(r'Photon energy$\: \rm [eV]$')
    plt.ylabel(r'Inelasticity')
    plt.savefig(f"{FIGURES_DIR}/inelasticity.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/inelasticity.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_photopion_inelasticity()

# ----------------------------------------------------------------------------------------------------