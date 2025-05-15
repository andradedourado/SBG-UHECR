from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'legend.fontsize': 'large',
'legend.title_fontsize': 'large',
'axes.labelsize': 'x-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

DATA_DIR = "../data"
FIGURES_DIR = "../figures"
RESULTS_DIR = "../results"

Auger_data = np.loadtxt(f"{DATA_DIR}/CombinedSpectrum_ICRC2019.txt")

cm_to_km = 1e-5
eV_to_erg = 1.602176634e-12

c = 2.99792458e10 # cm / s

# ----------------------------------------------------------------------------------------------------
def norm(E, J):

    f = interp1d(E, J, kind = 'cubic')
    norm = Auger_data[-6,1] * 0.1 / f(10**Auger_data[-6,0]) # eV / s
    
    np.savetxt(f"{RESULTS_DIR}/sbg_luminosity.dat", [norm], fmt = "%.15e")
    print(norm * eV_to_erg) # erg / s
    return norm    

# ----------------------------------------------------------------------------------------------------
def plot_Auger_intensity():

    threshold = 38e18 # eV
    mask = 10**Auger_data[:,0] > threshold

    E = 10**Auger_data[:,0]
    J = Auger_data[:,1]
    J_lower = Auger_data[:,2]
    J_upper = Auger_data[:,3]

    J_err = [J - J_lower, J_upper - J] 
    J_err_scaled = [E**2 * (J_err[0]), E**2 * (J_err[1])]
    plt.errorbar(np.log10(E), E**2 * J, yerr = J_err_scaled, color = 'k', linestyle = 'none', marker = '.', label = 'Auger (ICRC 2019)')
    plt.errorbar(np.log10(E[mask]), E[mask]**2 * J[mask] * 0.1, yerr = J_err_scaled[0][mask] * 0.1, color = 'gray', linestyle = 'none', marker = '.')

# ----------------------------------------------------------------------------------------------------
def plot_intensity():

    J = np.zeros(100)

    for igal in range(44):
        data = np.loadtxt(f"{RESULTS_DIR}/uhecr_density_{igal+2:02d}_injSBG.dat")
        J += c / (4 * np.pi) * data[:,1] / cm_to_km**2 # km^-2 eV^-2 sr^-1

    mask = J > 0
    plt.plot(np.log10(data[:,0][mask]), data[:,0][mask]**2 * J[mask] * norm(data[:,0], J) , color = 'k')
    
    plot_Auger_intensity()
    
    plt.yscale('log')
    plt.xlabel(r'$\log_{10}(\rm Energy/eV)$')
    plt.ylabel(r'$E^2 \times {\rm Intensity} \: \rm [eV \, km^{-2} \, sr^{-1} \, yr^{-1}]$')
    plt.savefig(f"{FIGURES_DIR}/normalized_intensity_starburst_galaxies.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/normalized_intensity_starburst_galaxies.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    plot_intensity()

# ----------------------------------------------------------------------------------------------------