from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'legend.fontsize': 'large',
'legend.title_fontsize': 'large',
'axes.labelsize': 'x-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

REFERENCES_DIR = "../references"
FIGURES_DIR = "../figures"
RESULTS_DIR = "../results"

# ----------------------------------------------------------------------------------------------------
def get_reference_figure_and_set_labels(l):

    if l == 'gmm':
        plt.xlabel(r'$x = E_\gamma / E_p$')
        plt.ylabel(r'$x^2 F_\gamma (x, E_p)$')
        return 'Figure6'
    
    elif l == 'e':
        plt.xlabel(r'$x = E_e / E_p$')
        plt.ylabel(r'$x^2 F_e (x, E_p)$')
        return 'Figure8'
    
    elif l == 'nu_mu':
        plt.xlabel(r'$x = E_{\nu_\mu} / E_p$')
        plt.ylabel(r'$x F_{\nu_\mu} (x, E_p)$')
        return 'Figure9'

# ----------------------------------------------------------------------------------------------------
def get_energy_string_and_column_index(Ep):

    if Ep == 1:
        return '1TeV', 1

    elif Ep == 30:
        return '30TeV', 2
    
    elif Ep == 300:
        return '300TeV', 3
    
    elif Ep == 3000:
        return '3000TeV', 4
    
    elif Ep == 0.1:
        return '0_1TeV', 1
    
    elif Ep == 100:
        return '100TeV', 2
    
    elif Ep == 1000:
        return '1000TeV', 3

# ----------------------------------------------------------------------------------------------------
def plot_KAB06_pp(l, Ep):

    data_KAB06 = np.loadtxt(f"{REFERENCES_DIR}/KAB06_{get_reference_figure_and_set_labels(l)}_{get_energy_string_and_column_index(Ep)[0]}.dat")
    data_LAD = np.loadtxt(f"{RESULTS_DIR}/KAB06_xchecks_F_{l}.dat")

    at = AnchoredText(r'$E_p = {0} \: \rm TeV$'.format(Ep), loc = 'upper left', frameon = False, prop = {'fontsize': 'x-large'})
    plt.gca().add_artist(at)

    plt.plot(data_KAB06[:,0], data_KAB06[:,1], color = 'k', ls = '--', label = 'KAB06')
    plt.plot(data_LAD[:,0], data_LAD[:,0]**2 * data_LAD[:,get_energy_string_and_column_index(Ep)[1]], color = 'r', ls = '-', label = 'LAD')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1e-3, 1e0])
    plt.ylim([1e-3, 1e-1])
    plt.legend(loc = 'upper right')
    plt.savefig(f"{FIGURES_DIR}/KAB06_F_{l}_comparison_{get_energy_string_and_column_index(Ep)[0]}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/KAB06_F_{l}_comparison_{get_energy_string_and_column_index(Ep)[0]}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # plot_KAB06_pp('gmm', 1)
    # plot_KAB06_pp('gmm', 30)
    # plot_KAB06_pp('gmm', 300)
    # plot_KAB06_pp('gmm', 3000)

    # plot_KAB06_pp('e', 0.1)
    # plot_KAB06_pp('e', 100)
    # plot_KAB06_pp('e', 1000)

    plot_KAB06_pp('nu_mu', 0.1)
    plot_KAB06_pp('nu_mu', 100)
    plot_KAB06_pp('nu_mu', 1000)

# ----------------------------------------------------------------------------------------------------
