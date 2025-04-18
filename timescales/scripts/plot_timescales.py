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

PARTICLES_LEGEND = [r'$^1\mathrm{H}$', r'$^4\mathrm{He}$', r'$^{14}\mathrm{N}$', r'$^{28}\mathrm{Si}$', r'$^{56}\mathrm{Fe}$']

# ----------------------------------------------------------------------------------------------------
def plot_timescales():

    data_advection = np.loadtxt(f"{RESULTS_DIR}/timescales_advection.dat") 
    data_diff_1H = np.loadtxt(f"{RESULTS_DIR}/timescales_diff_1H.dat")
    data_diff_56Fe = np.loadtxt(f"{RESULTS_DIR}/timescales_diff_56Fe.dat")
    data_spal_1H = np.loadtxt(f"{RESULTS_DIR}/timescales_spal_1H.dat")
    data_spal_56Fe = np.loadtxt(f"{RESULTS_DIR}/timescales_spal_56Fe.dat")

    plt.plot(np.log10(data_advection[:,0] * 1.e18), data_advection[:,1], label = 'Advection')
    plt.plot(np.log10(data_diff_1H[:,0] * 1.e18), data_diff_1H[:,1], label = f'Diffusion ({PARTICLES_LEGEND[0]})')
    plt.plot(np.log10(data_diff_56Fe[:,0] * 1.e18), data_diff_56Fe[:,1], label = f'Diffusion ({PARTICLES_LEGEND[4]})')
    plt.plot(np.log10(data_spal_1H[:,0] * 1.e18), data_spal_1H[:,1], label = f'Hadronic interaction ({PARTICLES_LEGEND[0]})')
    plt.plot(np.log10(data_spal_56Fe[:,0] * 1.e18), data_spal_56Fe[:,1], label = f'Hadronic interaction ({PARTICLES_LEGEND[4]})')
    plt.yscale('log')
    plt.ylim([1.e2, 1.e7])
    plt.xlabel(r'$\log_{10}{({\rm Energy}/{\rm eV})}$')
    plt.ylabel(r'Timescales$\: \rm [yr]$')
    plt.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 'lower left', ncols = 2, mode = "expand", borderaxespad = 0.)
    plt.savefig(f"{FIGURES_DIR}/timescales.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/timescales.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_timescales()

# ----------------------------------------------------------------------------------------------------