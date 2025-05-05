from matplotlib.offsetbox import AnchoredText
from matplotlib.ticker import FuncFormatter
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

gal_data = np.genfromtxt('starburst_galaxies.dat', dtype = None, encoding = None)

# ----------------------------------------------------------------------------------------------------
def get_spectrum_color(photon_spectrum):

    if photon_spectrum == 'IR':
        return '#FF5555'
    
    elif photon_spectrum == 'OPT':
        return '#5555FF'

# ----------------------------------------------------------------------------------------------------
def set_xlabel_and_get_index(param):

    if param == 'Cdil':
        plt.ticklabel_format(style = 'sci', scilimits = (0,0), axis = 'x', useMathText = True)
        plt.xlabel(r'Normalization')
        return 0

    elif param == 'temp':
        plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.0f}' if x.is_integer() else f'{x}'))
        plt.xlabel(r'Temperature$\: \rm [K]$')
        return 1

    elif param == 'sgm':
        plt.xlabel(r'$\sigma$')
        return 2

# ----------------------------------------------------------------------------------------------------
def plot_SED_fit_histogram(photon_spectrum, param):

    hist = []
    index = set_xlabel_and_get_index(param)

    for igal, galaxy in enumerate(gal_data):
        if galaxy[2] == 'NGC4038/9':
            data_NGC4038 = np.loadtxt(f"{RESULTS_DIR}/fitted_params_36_NGC4038_{photon_spectrum}.dat") 
            data_NGC4039 = np.loadtxt(f"{RESULTS_DIR}/fitted_params_36_NGC4039_{photon_spectrum}.dat")
            hist.append(data_NGC4038[index])
            hist.append(data_NGC4039[index])
        else:
            data = np.loadtxt(f"{RESULTS_DIR}/fitted_params_{(igal + 2):02d}_{galaxy[2]}_{photon_spectrum}.dat")
            hist.append(data[index])

    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.0f}' if y.is_integer() else f'{y}'))

    plt.gca().add_artist(AnchoredText(f'{photon_spectrum}', loc = 'upper center', frameon = False, prop = {'fontsize': 'x-large'}))

    plt.hist(hist, color = get_spectrum_color(photon_spectrum), edgecolor = 'black')
    plt.ylabel(r'Number of galaxies')
    plt.savefig(f"{FIGURES_DIR}/hist_{param}_{photon_spectrum}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/hist_{param}_{photon_spectrum}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_SED_fit_histogram('IR', 'temp')
    plot_SED_fit_histogram('OPT', 'temp')

    plot_SED_fit_histogram('IR', 'Cdil')
    plot_SED_fit_histogram('OPT', 'Cdil')

    plot_SED_fit_histogram('IR', 'sgm')

# ----------------------------------------------------------------------------------------------------