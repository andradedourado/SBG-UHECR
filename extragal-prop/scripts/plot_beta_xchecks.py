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

cm_to_Mpc = 1 / 3.086e24 
s_to_yr = 3.17e-8

c = 3e10 * cm_to_Mpc # Mpc / s

# ----------------------------------------------------------------------------------------------------
def plot_beta_xchecks():

    data_AC = np.loadtxt(f"{REFERENCES_DIR}/AC_log10_beta.dat")
    data_LAD = np.loadtxt(f"{RESULTS_DIR}/beta.dat")
    data_RA = np.loadtxt(f"{REFERENCES_DIR}/RA_attenuation_length_1H.dat")
    
    plt.plot(data_AC[:,0], pow(10, data_AC[:,1])**-1 * c / s_to_yr, label = 'AC')
    plt.plot(np.log10(data_LAD[:,0]), (data_LAD[:,1])**-1 * c, label = 'LAD')
    plt.plot(np.log10(data_RA[:,0]), data_RA[:,1], label = 'RA')

    plt.xlabel(r'$\log_{10}{({\rm Energy}/{\rm eV})}$')
    plt.ylabel(r'Attenuation length$\: \rm [Mpc]$')
    plt.yscale('log')
    plt.legend(title = 'Results')
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_beta_xchecks()

# ----------------------------------------------------------------------------------------------------