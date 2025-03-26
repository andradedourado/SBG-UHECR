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

# ----------------------------------------------------------------------------------------------------
def plot_transport_solution():

    data = np.loadtxt(f"{RESULTS_DIR}/transport_sol_1H.dat")

    plt.plot(data[:,0], data[:,1])
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r'Energy$\: \rm [eV]$')
    plt.ylabel(r'$n(E) \: \rm [cm^{-3}]$')
    plt.ylim([1e10, 1e14])
    plt.savefig(f"{FIGURES_DIR}/transport_sol_1H.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/transport_sol_1H.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_transport_solution()

# ----------------------------------------------------------------------------------------------------