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
 
c = 3.07e-7 # Mpc / yr
mp = 1e9 # eV

# ----------------------------------------------------------------------------------------------------
def plot_phopion_timescales_xchecks():

    data_CE = np.loadtxt(f"{REFERENCES_DIR}/CE_timescales_photopion_CMB.dat")
    data_LAD = np.loadtxt(f"{RESULTS_DIR}/timescales_photopion_1H_CMB.dat")

    at = AnchoredText(r'P$\rm \pi$ | CMB', loc = 'upper center', frameon = False, prop = {'fontsize': 'x-large'})
    plt.gca().add_artist(at)
    
    plt.plot(data_CE[:,0] * mp, data_CE[:,1] / c, label = 'CE')
    plt.plot(data_LAD[:,0], data_LAD[:,1], label = 'LAD')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Energy$\: \rm [eV]$')
    plt.ylabel(r'Timescales$\: \rm [yr]$')
    plt.legend(title = 'Results')
    plt.savefig(f"{FIGURES_DIR}/timescales_photopion_xchecks_CMB.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/timescales_photopion_xchecks_CMB.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
def plot_pp_timescales_xchecks():

    data_AA = np.loadtxt(f"{REFERENCES_DIR}/AA_timescales_pp.dat")
    data_AC_IR = np.loadtxt(f"{REFERENCES_DIR}/AC_timescales_pp_IR.txt")
    data_AC_OPT = np.loadtxt(f"{REFERENCES_DIR}/AC_timescales_pp_OPT.txt")
    data_LAD = np.loadtxt(f"{RESULTS_DIR}/timescales_pairproduction_1H.dat")

    at = AnchoredText('PP', loc = 'upper center', frameon = False, prop = {'fontsize': 'x-large'})
    plt.gca().add_artist(at)

    plt.plot(data_AA[:,0], data_AA[:,1], label = 'AA')
    plt.plot(data_AC_IR[:,0], (data_AC_IR[:,1] + data_AC_OPT[:,1])**-1, label = 'AC')
    plt.plot(data_LAD[:,0], data_LAD[:,1], label = 'LAD')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Energy$\: \rm [eV]$')
    plt.ylabel(r'Timescales$\: \rm [yr]$')
    plt.legend(title = "Results")
    plt.savefig(f"{FIGURES_DIR}/timescales_pairproducton_xchecks.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/timescales_pairproducton_xchecks.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_pp_timescales_xchecks()
    plot_phopion_timescales_xchecks()

# ----------------------------------------------------------------------------------------------------