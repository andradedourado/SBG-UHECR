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
def plot_timescales_xchecks(inelasticity):

    data_photopion_LAD = np.loadtxt(f"{RESULTS_DIR}/timescales_photopion_1H.dat")
    data_pp_LAD = np.loadtxt(f"{RESULTS_DIR}/timescales_pairproduction_1H.dat")
    data_LAD = (1 / data_photopion_LAD[:,1] + 1 / data_pp_LAD[:,1])**-1
    
    data_Condo = np.loadtxt(f"{REFERENCES_DIR}/timescales_Condorelli_1H.dat")

    if inelasticity == True:
        data_photopion_Condo_LAD = np.loadtxt(f"{RESULTS_DIR}/timescales_photopion_xchecks_1H_wInelasticity.dat")
    elif inelasticity == False:
        data_photopion_Condo_LAD = np.loadtxt(f"{RESULTS_DIR}/timescales_photopion_xchecks_1H_woInelasticity.dat")
    data_pp_Condo_LAD = np.loadtxt(f"{RESULTS_DIR}/timescales_pp_xchecks_1H.dat")
    data_Condo_LAD = (1 / data_photopion_Condo_LAD[:,1] + 1 / data_pp_Condo_LAD[:,1])**-1

    plt.plot(data_Condo[:,0], data_Condo[:,1], label = 'Condo+23')
    plt.plot(np.log10(data_photopion_Condo_LAD[:,0]), data_Condo_LAD, label = 'Condo+23 by LAD')
    plt.plot(np.log10(data_photopion_LAD[:,0]), data_LAD, label = 'LAD')
    plt.yscale('log')
    plt.xlabel(r'$\log_{10}{({\rm Energy}/{\rm eV})}$')
    plt.ylabel(r'Timescales$\: \rm [yr]$')
    plt.legend(title = "Results")
    
    if inelasticity == True:
        plt.savefig(f"{FIGURES_DIR}/timescales_xchecks_wInelasticity.pdf", bbox_inches = 'tight')
        plt.savefig(f"{FIGURES_DIR}/timescales_xchecks_wInelasticity.png", bbox_inches = 'tight', dpi = 300)
    elif inelasticity == False:
        plt.savefig(f"{FIGURES_DIR}/timescales_xchecks_woInelasticity.pdf", bbox_inches = 'tight')
        plt.savefig(f"{FIGURES_DIR}/timescales_xchecks_woInelasticity.png", bbox_inches = 'tight', dpi = 300)
    
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_timescales_xchecks(True)
    plot_timescales_xchecks(False)

# ----------------------------------------------------------------------------------------------------