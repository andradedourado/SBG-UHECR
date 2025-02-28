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

PARTICLES_LEGEND = [r'$^{14}\mathrm{N}$', r'$^{28}\mathrm{Si}$', r'$^{56}\mathrm{Fe}$']
PARTICLES = ['14N', '28Si', '56Fe']
ZS = [7, 14, 26]

# ----------------------------------------------------------------------------------------------------
def iZ(Z):

    try:
        return ZS.index(Z)
    except ValueError:
        raise ValueError(f"Z ({Z}) not found in ZS.")

# ----------------------------------------------------------------------------------------------------
def get_color(Z):

    if Z == 7:
        return 'g'
    elif Z == 14:
        return 'r'
    elif Z == 26:
        return 'orange'
    else: 
        raise ValueError(f"Unknown Z value: {Z}") 

# ----------------------------------------------------------------------------------------------------
def plot_timescales():

    data_photopion_1H = np.loadtxt(f"{RESULTS_DIR}/timescales_photopion_1H.dat")
    data_pp_1H = np.loadtxt(f"{RESULTS_DIR}/timescales_pairproduction_1H.dat")    
    plt.plot(np.log10(data_photopion_1H[:,0]), data_photopion_1H[:,1], ls = '-.', color = 'b')
    plt.plot(np.log10(data_pp_1H[:,0]), data_pp_1H[:,1], ls = ':', color = 'b')
    plt.plot(np.log10(data_photopion_1H[:,0]), (1 / data_photopion_1H[:,1] + 1 / data_pp_1H[:,1])**-1, ls = '-', color = 'b')

    for Z in ZS:

        data_pd = np.loadtxt(f"{RESULTS_DIR}/timescales_photodisintegration_{PARTICLES[iZ(Z)]}.dat")
        data_photopion = np.loadtxt(f"{RESULTS_DIR}/timescales_photopion_{PARTICLES[iZ(Z)]}.dat")
        data_pp = np.loadtxt(f"{RESULTS_DIR}/timescales_pairproduction_{PARTICLES[iZ(Z)]}.dat")

        color = get_color(Z)

        plt.plot(np.log10(data_pd[:,0]), data_pd[:,1], ls = '--', color = color)
        plt.plot(np.log10(data_photopion[:,0]), data_photopion[:,1], ls = '-.', color = color)
        plt.plot(np.log10(data_pp[:,0]), data_pp[:,1], ls = ':', color = color)
        plt.plot(np.log10(data_pd[:,0]), (1 / data_pd[:,1] + 1 / data_photopion[:,1] + 1 / data_pp[:,1])**-1, ls = '-', color = color)

    H_label = lines.Line2D([], [], color = 'blue', label = r'$^{1}$H')
    N_label = lines.Line2D([], [], color = 'green', label = r'$^{14}$N')
    Si_label = lines.Line2D([], [], color = 'red', label = r'$^{28}$Si')
    Fe_label = lines.Line2D([], [], color = 'orange', label = r'$^{56}$Fe')
    lgnd = plt.legend(title = 'Nucleus', handles = [H_label, N_label, Si_label, Fe_label], frameon = True, loc = 'upper right')
    plt.gca().add_artist(lgnd)

    photodisintegration_label = lines.Line2D([], [], color = 'black', ls = '--', label = r'$\rm PD$')
    photopion_label = lines.Line2D([], [], color = 'black', ls = '-.', label = r'$\rm P\pi$')
    pairproduction_label = lines.Line2D([], [], color = 'black', ls = ':', label = r'$\rm PP$')
    total_label = lines.Line2D([], [], color = 'black', ls = '-', label = 'Total')
    lgnd = plt.legend(title = 'Interaction', handles = [photodisintegration_label, photopion_label, pairproduction_label, total_label], frameon = True, loc = 'lower left')
    plt.gca().add_artist(lgnd)

    plt.yscale('log')
    plt.ylim(top = 1e8)
    plt.xlabel(r'$\log_{10}{({\rm Energy}/{\rm eV})}$')
    plt.ylabel(r'Timescales$\: \rm [yr]$')
    plt.savefig(f"{FIGURES_DIR}/photohadronic_timescales.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/photohadronic_timescales.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_timescales()

# ----------------------------------------------------------------------------------------------------