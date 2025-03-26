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
def get_figure_and_set_labels(l):

    if l == 'gmm':
        plt.xlabel(r'$x = E_\gamma / E_p$')
        plt.ylabel(r'$x \Phi_\gamma (\eta, x)$')
        return 'Figure2'
    
    elif l == 'pos':
        plt.xlabel(r'$x = E_{e^+} / E_p$')
        plt.ylabel(r'$x \Phi_{e^+} (\eta, x)$')
        return 'Figure3'
    
    elif l == 'anu_mu':
        plt.xlabel(r'$x = E_{\bar{\nu}_\mu} / E_p$')
        plt.ylabel(r'$x \Phi_{\bar{\nu}_\mu} (\eta, x)$')
        return 'Figure4'
    
    elif l == 'nu_mu':
        plt.xlabel(r'$x = E_{\nu_\mu} / E_p$')
        plt.ylabel(r'$x \Phi_{\nu_\mu} (\eta, x)$')
        return 'Figure5'

# ----------------------------------------------------------------------------------------------------
def set_limit_and_get_position_and_index(l, eta_over_eta0):

    if eta_over_eta0 == 1.5:
        plt.ylim([1.e-21, 1.e-18])
        return 'Left', 5
    
    elif eta_over_eta0 == 30:
        plt.ylim([1.e-20, 1.e-17])
        if l == 'gmm':
            return 'Right', 20
        else:
            return 'Right', 19

# ----------------------------------------------------------------------------------------------------
def plot_KA08_parametrizations(l, eta_over_eta0):

    data_KA08 = np.loadtxt(f"{REFERENCES_DIR}/KA08_{get_figure_and_set_labels(l)}_{set_limit_and_get_position_and_index(l, eta_over_eta0)[0]}.dat")
    data_LAD = np.loadtxt(f"{RESULTS_DIR}/KA08_phi_{l}.dat")

    at = AnchoredText(r'$\eta = {0} \eta_0$'.format(eta_over_eta0), loc = 'upper left', frameon = False, prop = {'fontsize': 'x-large'})
    plt.gca().add_artist(at)

    plt.plot(data_KA08[:,0], data_KA08[:,1], color = 'k', ls = '--', label = 'KA08')
    plt.plot(data_LAD[:,0], data_LAD[:,0] * data_LAD[:,set_limit_and_get_position_and_index(l, eta_over_eta0)[1]], color = 'r', ls = '-', label = 'LAD')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.savefig(f"{FIGURES_DIR}/KA08_phi_{l}_comparison_{set_limit_and_get_position_and_index(l, eta_over_eta0)[1]:02d}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/KA08_phi_{l}_comparison_{set_limit_and_get_position_and_index(l, eta_over_eta0)[1]:02d}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_KA08_parametrizations('gmm', 1.5)
    plot_KA08_parametrizations('gmm', 30)
    plot_KA08_parametrizations('pos', 1.5)
    plot_KA08_parametrizations('pos', 30)
    plot_KA08_parametrizations('anu_mu', 1.5)
    plot_KA08_parametrizations('anu_mu', 30)
    plot_KA08_parametrizations('nu_mu', 1.5)
    plot_KA08_parametrizations('nu_mu', 30)

# ----------------------------------------------------------------------------------------------------
