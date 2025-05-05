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

# ----------------------------------------------------------------------------------------------------
def get_zg_str(zg):

    if zg == int(zg):
        return f"{int(zg)}"
    else:
        return str(zg).replace(".", "_")

# ----------------------------------------------------------------------------------------------------
def plot_single_source_solution():

    for izg, zg in enumerate([0.05, 0.5, 1, 2]):
        data = np.loadtxt(f"{RESULTS_DIR}/single_source_solution_zg{get_zg_str(zg)}.dat")
        mask = (data[:,1] < 1e30) & np.isfinite(data[:,1])
        plt.plot(data[:,0][mask], data[:,1][mask], label = f'${zg}$')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r"Energy$\: \rm [eV]$")
    plt.ylabel(r"$Q_0(E_{\rm g}(E, z_{\rm g}))dE_{\rm g}(E, z_{\rm g})/dE \: \rm [arb. units]$")
    plt.legend(title = r'$z_{\rm g}$')
    plt.savefig(f"{FIGURES_DIR}/single_source_solution.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/single_source_solution.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_single_source_solution()

# ----------------------------------------------------------------------------------------------------