from matplotlib.offsetbox import AnchoredText
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
def get_zg_arr(zg):

    if zg == int(zg):
        return f"{int(zg)}"
    else:
        return str(zg).replace(".", "_")
    
# ----------------------------------------------------------------------------------------------------
def plot_dEg_dE_vs_E(izg, zg):

    data_AC = np.loadtxt(f"{REFERENCES_DIR}/AC_dEgdE_vs_E.txt")
    data_LAD = np.loadtxt(f"{RESULTS_DIR}/dEg_dE_zg{get_zg_arr(zg)}.dat")

    plt.plot(data_AC[:,0], data_AC[:,izg+1], label = 'AC')
    plt.plot(data_LAD[:,0], data_LAD[:,1], label = 'LAD')

    plt.gca().add_artist(AnchoredText(r'$z_{{\rm g}} = {}$'.format(zg), loc = 'lower right', frameon = False, prop = {'fontsize': 'x-large'}))

    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(top = 1e30)
    plt.xlabel(r"Energy$\: \rm [eV]$")
    plt.ylabel(r"$dE_{\rm g} (E, z_g) / dE$")
    plt.legend(title = 'Results')
    plt.savefig(f"{FIGURES_DIR}/dEg_dE_vs_E_xchecks_zg{get_zg_arr(zg)}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/dEg_dE_vs_E_xchecks_zg{get_zg_arr(zg)}.png", bbox_inches = 'tight', dpi = 300)
    plt.show() 

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for izg, zg in enumerate([0.05, 0.5, 1]):
        plot_dEg_dE_vs_E(izg, zg)

# ----------------------------------------------------------------------------------------------------