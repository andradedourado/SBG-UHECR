from matplotlib import lines
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
def plot_diffusion_timescales():

    data_MH19_diff_1H = np.loadtxt(f"{RESULTS_DIR}/timescales_MR19_diff_1H.dat")
    data_MH19_diff_56Fe = np.loadtxt(f"{RESULTS_DIR}/timescales_MR19_diff_56Fe.dat")
    data_diff_1H = np.loadtxt(f"{RESULTS_DIR}/timescales_diff_1H.dat")
    data_diff_56Fe = np.loadtxt(f"{RESULTS_DIR}/timescales_diff_56Fe.dat")
    
    plt.plot(data_MH19_diff_1H[:,0] * 1e18, data_MH19_diff_1H[:,1], color = 'b', ls = '--')
    plt.plot(data_MH19_diff_56Fe[:,0] * 1e18, data_MH19_diff_56Fe[:,1], color = 'orange', ls = '--')
    plt.plot(data_diff_1H[:,0] * 1e18, data_diff_1H[:,1], color = 'b', ls = '-.')
    plt.plot(data_diff_56Fe[:,0] * 1e18, data_diff_56Fe[:,1], color = 'orange', ls = '-.')

    H_label = lines.Line2D([], [], color = 'blue', label = r'$^{1}$H')
    Fe_label = lines.Line2D([], [], color = 'orange', label = r'$^{56}$Fe')
    lgnd = plt.legend(title = 'Nucleus', handles = [H_label, Fe_label], frameon = True, loc = 'upper right')
    plt.gca().add_artist(lgnd)

    Condo_label = lines.Line2D([], [], color = 'k', ls = '--', label = 'Condo+23')
    MH19_label = lines.Line2D([], [], color = 'k', ls = '-.', label = 'MR19')
    lgnd = plt.legend(title = r'$D(E)$', handles = [Condo_label, MH19_label], frameon = True, loc = 'upper right', bbox_to_anchor = (0.79, 1.))
    plt.gca().add_artist(lgnd)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Energy$\: \rm [eV]$')
    plt.ylabel(r'Timescales$\: \rm [yr]$')
    plt.savefig(f"{FIGURES_DIR}/diffusion_timescales.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/diffusion_timescales.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_diffusion_timescales()

# ----------------------------------------------------------------------------------------------------