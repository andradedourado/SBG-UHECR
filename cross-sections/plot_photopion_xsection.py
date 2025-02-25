import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({'legend.fontsize': 'large',
'legend.title_fontsize': 'large',
'axes.labelsize': 'x-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

mp = 0.938 # GeV 

# ----------------------------------------------------------------------------------------------------
def plot_photopion_xsections():

    s = []
    cross_section = []

    with open("xsecs_photopion_proton_sophia.txt", 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            
            values = line.strip().split(',')
            s.append(float(values[0]))  
            cross_section.append(float(values[1]))

    s = np.array(s)
    cross_section = np.array(cross_section)
    eps = (s - mp**2) / (2*mp) # s = m_p^2 + 2 m_p \epsilon'

    plt.plot(eps, cross_section, color = 'blue', linestyle = '-', label = 'SOPHIA')
    plt.xscale('log')
    plt.xlabel(r'Photon energy$\: \rm [GeV]$')
    plt.ylabel(r'Total inelastic cross section$\: \rm [mb]$')
    plt.legend()
    plt.savefig(f"photopion_xsection.pdf", bbox_inches = 'tight')
    plt.savefig(f"photopion_xsection.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_photopion_xsections()

# ----------------------------------------------------------------------------------------------------
