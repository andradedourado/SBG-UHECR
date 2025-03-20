from matplotlib import lines
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

ls = ['gmm', 'nu_mu', 'e']

# ----------------------------------------------------------------------------------------------------
def get_label(l):

    if l == 'gmm':
        return r'$\gamma$' 
    
    elif l == 'nu_mu':
        return r'$\nu_\mu$'
    
    elif l == 'e':
        return r'$e$'

# ----------------------------------------------------------------------------------------------------
def get_filename_parts_and_set_ylim(alp):

    if alp == 2:
        plt.ylim([1e-19, 4e-17])
        return 'Left', '2'

    elif alp == 1.5:
        plt.ylim([3e-19, 1e-16])
        return 'Right', '1_5'

# ----------------------------------------------------------------------------------------------------
def get_color(l):

    if l == 'gmm':
        return 'b'

    elif l == 'nu_mu':
        return 'g'

    elif l == 'e':
        return 'r'

# ----------------------------------------------------------------------------------------------------
def plot_KAB06_spectrum(alp):

    for l in ls:

        data_KAB06 = np.loadtxt(f"{REFERENCES_DIR}/KAB06_Figure12_{get_filename_parts_and_set_ylim(alp)[0]}_{l}.dat")
        data_LAD = np.loadtxt(f"{RESULTS_DIR}/KAB06_spectrum_{get_filename_parts_and_set_ylim(alp)[1]}_{l}.dat")

        plt.plot(data_KAB06[:,0], data_KAB06[:,1], color = get_color(l), ls = '--')
        plt.plot(data_LAD[:,0], data_LAD[:,0]**2 * data_LAD[:,1], color = get_color(l), label = get_label(l))

    at = AnchoredText(r'$\alpha = {0}$'.format(alp), loc = 'upper right', frameon = False, prop = {'fontsize': 'x-large'})
    plt.gca().add_artist(at)

    KAB06_label = lines.Line2D([], [], color = 'k', ls = '--', label = 'KAB06')
    LAD_label = lines.Line2D([], [], color = 'k', label = 'LAD')
    lgnd = plt.legend(title = 'Results', handles = [KAB06_label, LAD_label], frameon = True, loc = 'upper left')
    plt.gca().add_artist(lgnd)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r"Energy$\: \rm [TeV]$")
    plt.ylabel(r"$E^2 dN/dE \: \rm [TeV^{-1} \, cm^{-3} \, s^{-1}]$")
    plt.legend(title = 'Particle', loc = 'upper left', bbox_to_anchor = (0., 0.766))
    plt.savefig(f"{FIGURES_DIR}/KAB06_spectrum_{get_filename_parts_and_set_ylim(alp)[1]}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/KAB06_spectrum_{get_filename_parts_and_set_ylim(alp)[1]}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_KAB06_spectrum(2)
    plot_KAB06_spectrum(1.5)

# ----------------------------------------------------------------------------------------------------