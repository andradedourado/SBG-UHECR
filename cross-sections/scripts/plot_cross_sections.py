from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import ast
import numpy as np
import subprocess

plt.rcParams.update({'legend.fontsize': 'large',
'legend.title_fontsize': 'large',
'axes.labelsize': 'x-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

FIGURES_DIR = "../figures"

PARTICLES = ['1H', '4He', '14N', '28Si', '56Fe']
PARTICLES_LEGEND = [r'$^1\mathrm{H}$', r'$^4\mathrm{He}$', r'$^{14}\mathrm{N}$', r'$^{28}\mathrm{Si}$', r'$^{56}\mathrm{Fe}$']
ZS = [1, 2, 7, 14, 26]

GeV_to_MeV = 1e3
mp = 0.938 # GeV 

# ----------------------------------------------------------------------------------------------------
def iZ(Z):

    try:
        return ZS.index(Z)
    except ValueError:
        raise ValueError(f"Z ({Z}) not found in ZS.")

# ----------------------------------------------------------------------------------------------------
def execute_get_cross_section_TENDL2023(A, Z):

    output = subprocess.run(
        ['python3', f"get_cross_section_TENDL-2023.py", str(A), str(Z)],
        capture_output = True,
        text = True
    )
    
    if output.returncode != 0:
            raise RuntimeError(f'Error executing script: {output.stderr.strip()}')
    
    eps, cross_section = [ast.literal_eval(line) for line in output.stdout.strip().split('\n')]
    
    eps, cross_section = np.array(eps), np.array(cross_section)
    mask = cross_section > 0
    eps = eps[mask]
    cross_section = cross_section[mask] 

    return eps, cross_section

# ----------------------------------------------------------------------------------------------------
def get_photopion_xsections():

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

    return eps, cross_section 

# ----------------------------------------------------------------------------------------------------
def plot_xsections(A, Z):

    eps_pd, cross_section_pd = execute_get_cross_section_TENDL2023(A, Z)
    eps_photopion, cross_section_photopion = get_photopion_xsections()

    plt.plot(eps_pd, cross_section_pd, color = 'red', linestyle = '-', label = r'$\rm PD$')
    plt.plot(eps_photopion * GeV_to_MeV, cross_section_photopion * A, color = 'blue', linestyle = '-', label = r'$\rm P \pi$')

    at = AnchoredText(f'{PARTICLES_LEGEND[iZ(Z)]}', loc = 'upper left', frameon = False, prop = {'fontsize': 'x-large'})
    plt.gca().add_artist(at)
    
    plt.xscale('log')
    plt.xlabel(r'Photon energy$\: \rm [MeV]$')
    plt.ylabel(r'Total inelastic cross section$\: \rm [mb]$')
    plt.legend(title = 'Interactions')
    plt.savefig(f"{FIGURES_DIR}/xsections_{PARTICLES[iZ(Z)]}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/xsections_{PARTICLES[iZ(Z)]}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
def plot_photopion_xsections():

    eps, cross_section = get_photopion_xsections()

    plt.plot(eps, cross_section, color = 'blue', linestyle = '-', label = 'SOPHIA')
    plt.xscale('log')
    plt.xlabel(r'Photon energy$\: \rm [GeV]$')
    plt.ylabel(r'Total inelastic cross section$\: \rm [mb]$')
    plt.legend()
    plt.savefig(f"{FIGURES_DIR}/photopion_xsection.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/photopion_xsection.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_xsections(14, 7)
    plot_xsections(28, 14)
    plot_xsections(56, 26)
    # plot_photopion_xsections()

# ----------------------------------------------------------------------------------------------------
