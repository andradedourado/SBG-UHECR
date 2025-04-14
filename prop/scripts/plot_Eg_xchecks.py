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
def get_zg_str(zg):

    if zg == int(zg):
        return f"{int(zg)}"
    else:
        return str(zg).replace(".", "_")

# ----------------------------------------------------------------------------------------------------
def plot_Eg_vs_E_xchecks(izg, zg):

    data_AA = np.loadtxt(f"{REFERENCES_DIR}/AA_Eg_vs_E_zg{get_zg_str(zg)}.dat")
    data_AC = np.loadtxt(f"{REFERENCES_DIR}/AC_Eg_vs_E.txt")
    data_LAD = np.loadtxt(f"{RESULTS_DIR}/Eg_vs_E_xchecks_zg{get_zg_str(zg)}.dat")

    plt.plot(data_AA[:,0], data_AA[:,1], label = 'AA') 
    plt.plot(data_AC[:,0], data_AC[:,izg+1], label = 'AC')
    plt.plot(data_LAD[:,0], data_LAD[:,1], label = 'LAD') 

    plt.gca().add_artist(AnchoredText(r'$z_{{\rm g}} = {}$'.format(zg), loc = 'lower right', frameon = False, prop = {'fontsize': 'x-large'}))
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1e17, 1e22])
    plt.xlabel(r"Energy$\: \rm [eV]$")
    plt.ylabel(r"$E_{\rm g} \: \rm [eV]$")
    plt.legend(title = 'Results')
    plt.savefig(f"{FIGURES_DIR}/Eg_vs_E_xchecks_zg{get_zg_str(zg)}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/Eg_vs_E_xchecks_zg{get_zg_str(zg)}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
def plot_Eg_vs_E():

    for zg in [0.05, 0.5, 1, 2, 3]:
        data_LAD = np.loadtxt(f"{RESULTS_DIR}/Eg_vs_E_xchecks_zg{get_zg_str(zg)}.dat")
        plt.plot(data_LAD[:,0], data_LAD[:,1], label = f'{zg}') 

    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1e17, 1e22])
    plt.xlabel(r"Energy$\: \rm [eV]$")
    plt.ylabel(r"$E_{\rm g} \: \rm [eV]$")
    plt.legend(title = r'$z_{g}$')
    plt.savefig(f"{FIGURES_DIR}/Eg_vs_E.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/Eg_vs_E.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
def plot_Eg_vs_z_xchecks(iE, E):

    data_AC = np.loadtxt(f"{REFERENCES_DIR}/AC_Eg_vs_z.txt")
    data_LAD = np.loadtxt(f"{RESULTS_DIR}/Eg_vs_zg_xchecks_E{int(np.log10(E))}.dat")
    
    mask = (data_AC[:, iE+1] < 1e30) & np.isfinite(data_AC[:, iE+1])

    plt.plot(data_AC[:,0][mask], data_AC[:,iE+1][mask], label = 'AC')
    plt.plot(data_LAD[:,0], data_LAD[:,1], label = 'LAD')

    plt.gca().add_artist(AnchoredText(r'$E = 10^{{{}}} \: \rm eV$'.format(int(np.log10(E))), loc = 'upper center', frameon = False, prop = {'fontsize': 'x-large'}))

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r"Redshift")
    plt.ylabel(r"$E_{\rm g} \: \rm [eV]$")
    plt.legend(title = 'Results')
    plt.savefig(f"{FIGURES_DIR}/Eg_vs_z_xchecks_E{int(np.log10(E))}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/Eg_vs_z_xchecks_E{int(np.log10(E))}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for izg, zg in enumerate([0.05, 0.5]):
        plot_Eg_vs_E_xchecks(izg, zg)

    plot_Eg_vs_E()

    for iE, E in enumerate([1e17, 1e18, 1e19, 1e20, 1e21]):
        plot_Eg_vs_z_xchecks(iE, E)

# ----------------------------------------------------------------------------------------------------