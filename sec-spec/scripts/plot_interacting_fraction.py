from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np 

plt.rcParams.update({'legend.fontsize': 'large',
'legend.title_fontsize': 'large',
'axes.labelsize': 'x-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

FIGURES_DIR = "../figures"
TIMESCALES_DIR = "../../timescales/results"

yr_to_s = 60 * 60 * 24 * 365.25

# ----------------------------------------------------------------------------------------------------
def get_total_timescale(E):

    files = ["timescales_advection.dat", "timescales_diff_1H.dat", "timescales_pairproduction_1H.dat", "timescales_photopion_1H.dat", "timescales_spal_1H.dat"]
    
    energy_conversion = {"timescales_advection.dat": 1.e18, "timescales_diff_1H.dat": 1.e18, "timescales_pairproduction_1H.dat": 1.0, "timescales_photopion_1H.dat": 1.0, "timescales_spal_1H.dat": 1.e18}

    interpolated_tau = []

    for file in files:
        data = np.loadtxt(f"{TIMESCALES_DIR}/{file}")
        E_table = data[:,0] * energy_conversion[file]
        tau_table = data[:,1] * yr_to_s  

        f_interp = interp1d(E_table, tau_table, kind = 'cubic', bounds_error = False, fill_value = np.inf)
        interpolated_tau.append(f_interp(E))

    total_tau_inv = sum(1 / tau for tau in interpolated_tau)
    return np.where(total_tau_inv > 0, 1 / total_tau_inv, np.inf) 

# ----------------------------------------------------------------------------------------------------
def compute_interaction_timescales(E, interaction):

    if interaction == 'pp':
        data = np.loadtxt(f"{TIMESCALES_DIR}/timescales_spal_1H.dat")
        E_table = data[:,0] * 1.e18
        tau_table = data[:,1] * yr_to_s

        f_interp = interp1d(E_table, tau_table, kind = 'cubic', bounds_error = False, fill_value = np.inf)
        return f_interp(E)

    elif interaction == 'pgamma':
        data = np.loadtxt(f"{TIMESCALES_DIR}/timescales_photopion_1H.dat")
        E_table = data[:,0]
        tau_table = data[:,1] * yr_to_s

        f_interp = interp1d(E_table, tau_table, kind = 'cubic', bounds_error = False, fill_value = np.inf)
        return f_interp(E)

# ----------------------------------------------------------------------------------------------------
def plot_interacting_fraction():

    E = np.logspace(17, 21, num = 500)

    plt.plot(E, get_total_timescale(E) / compute_interaction_timescales(E, 'pp'), label = r'$pp$')
    plt.plot(E[:-4], (get_total_timescale(E) / compute_interaction_timescales(E, 'pgamma'))[:-4], label = r'$p\gamma$') 
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Energy$\: \rm [eV]$')
    plt.ylabel(r'$\tau_{\rm total} / \tau_{\rm interaction}$')
    plt.legend(title = 'Interaction')
    plt.savefig(f"{FIGURES_DIR}/interacting_fraction.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/interacting_fraction.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_interacting_fraction()

# ----------------------------------------------------------------------------------------------------