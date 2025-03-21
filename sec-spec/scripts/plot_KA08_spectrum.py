from matplotlib.offsetbox import AnchoredText
from matplotlib import lines
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

ls_Left = ['gmm', 'pos', 'e']
ls_Right = ['anu_mu', 'nu_mu', 'nu_e', 'anu_e']

# ----------------------------------------------------------------------------------------------------
def get_figure_and_set_labels(Ep):

    if Ep == 1e20:
        plt.ylim([1e-19, 2e-16])
        return 'Figure6'
    
    elif Ep == 1e21:
        plt.ylim([1e-18, 1e-14])
        return 'Figure7'

# ----------------------------------------------------------------------------------------------------
def get_color_and_label(l):

    if l == 'gmm':
        return 'tab:blue', r'$\gamma$'
    
    elif l == 'pos':
        return 'tab:orange', r'$e^+$'
    
    elif l == 'e':
        return 'tab:green', r'$e^-$'
    
    elif l == 'anu_mu':
        return 'tab:red', r'$\bar{\nu}_\mu$'
    
    elif l == 'nu_mu':
        return 'tab:purple', r'$\nu_\mu$'
    
    elif l == 'nu_e':
        return 'tab:brown', r'$\nu_e$'
    
    elif l == 'anu_e':
        return 'tab:pink', r'$\bar{\nu}_e$'

# ----------------------------------------------------------------------------------------------------
def plot_mono_spectrum(ls, Ep):

    for l in ls:
        data_KA08 = np.loadtxt(f"{REFERENCES_DIR}/KA08_{get_figure_and_set_labels(Ep)}_{l}.dat")
        data_LAD = np.loadtxt(f"{RESULTS_DIR}/KA08_mono_spectrum_CMB_Ep{int(np.log10(Ep)):02d}_{l}.dat")
        plt.plot(data_KA08[:,0], data_KA08[:,1], color = get_color_and_label(l)[0], ls = '--')
        plt.plot(data_LAD[:,0], data_LAD[:,0] * data_LAD[:,1], color = get_color_and_label(l)[0], label = f'{get_color_and_label(l)[1]}')

    at = AnchoredText(r'$E_p = 10^{{{0}}} \: \rm eV$'.format(int(np.log10(Ep))), loc = 'upper left', frameon = False, prop = {'fontsize': 'x-large'})
    plt.gca().add_artist(at)

    KA08_lgnd = lines.Line2D([], [], color = 'black', ls = '--', label = 'KA08')
    LAD_lgnd = lines.Line2D([], [], color = 'black', ls = '-', label = 'LAD')
    lgnd = plt.legend(title = 'Result', handles = [KA08_lgnd, LAD_lgnd], frameon = True, loc = 'upper right')
    plt.gca().add_artist(lgnd)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$x = E/E_p$')
    plt.ylabel(r'$x \, dN/dx \: \rm [s^{-1}]$')
    plt.legend(title = 'Particle', loc = 'upper right', bbox_to_anchor = (1., 0.767))

    if ls == ls_Left:
        plt.savefig(f"{FIGURES_DIR}/mono_spectrum_Ep{int(np.log10(Ep)):02d}_Left.pdf", bbox_inches = 'tight')
        plt.savefig(f"{FIGURES_DIR}/mono_spectrum_Ep{int(np.log10(Ep)):02d}_Left.png", bbox_inches = 'tight', dpi = 300)
    
    elif ls == ls_Right:
        plt.savefig(f"{FIGURES_DIR}/mono_spectrum_Ep{int(np.log10(Ep)):02d}_Right.pdf", bbox_inches = 'tight')
        plt.savefig(f"{FIGURES_DIR}/mono_spectrum_Ep{int(np.log10(Ep)):02d}_Right.png", bbox_inches = 'tight', dpi = 300)
    
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_mono_spectrum(ls_Left, 1e20)
    plot_mono_spectrum(ls_Right, 1e20)
    plot_mono_spectrum(ls_Left, 1e21)
    plot_mono_spectrum(ls_Right, 1e21)

# ----------------------------------------------------------------------------------------------------