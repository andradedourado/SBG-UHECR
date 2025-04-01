from scipy.integrate import simps
import ast
import numpy as np 
import subprocess

RESULTS_DIR = "../results"
XSECTIONS_DIR = "../../cross-sections/scripts"
PARTICLES = ['1H', '4He', '14N', '28Si', '56Fe']
ZS = [1, 2, 7, 14, 26]

mbarn_to_m2 = 1.e-31
GeV_to_eV = 1.e9
MeV_to_eV = 1.e6     
s_to_yr = 1 / (60 * 60 * 24 * 365.25) 

c =  299792458       # m/s
hbar = 6.5821220e-16 # eV.s
kB = 8.6173303e-5    # eV/K
mp = 0.9383e9        # eV 

T_IR = 3.5e-3 / kB    # K
T_OPT = 332.5e-3 / kB # K
NORM_IR = 0.15268288372409347
NORM_OPT = 2.810860522825165e-09

# ----------------------------------------------------------------------------------------------------
def iZ(Z):

    try:
        return ZS.index(Z)
    except ValueError:
        raise ValueError(f"Z ({Z}) not found in ZS.")

# ----------------------------------------------------------------------------------------------------
def execute_get_cross_section_TENDL2023(A, Z):

    output = subprocess.run(
        ['python3', f"{XSECTIONS_DIR}/get_cross_section_TENDL-2023.py", str(A), str(Z)],
        capture_output = True,
        text = True
    )
    
    if output.returncode != 0:
            raise RuntimeError(f'Error executing script: {output.stderr.strip()}')
    
    eps, cross_section = [ast.literal_eval(line) for line in output.stdout.strip().split('\n')]

    return np.array(eps), np.array(cross_section)

# ----------------------------------------------------------------------------------------------------
def I(eps, Gmm):

    return -kB / (np.pi**2 * (hbar * c)**3) * (T_IR * NORM_IR * np.log(1. - np.exp(-(eps)/(2 * Gmm * kB * T_IR))) + T_OPT * NORM_OPT * np.log(1. - np.exp(-(eps)/(2 * Gmm * kB * T_OPT))))

# ----------------------------------------------------------------------------------------------------
def inelasticity(x): 

    Y_inf = 0.47
    xb = 6e9 # eV
    dlt = 0.33
    s = 0.15

    return Y_inf * (x / xb)**dlt / (1 + (x / xb)**(dlt / s))**s

# ----------------------------------------------------------------------------------------------------
def get_eps_and_cross_sections(A, Z, interaction):

    if interaction == 'photodisintegration':
        eps, cross_section = execute_get_cross_section_TENDL2023(A, Z)
        mask = cross_section > 0
        eps = eps[mask] * MeV_to_eV 
        cross_section = cross_section[mask] * mbarn_to_m2
        return eps, cross_section
    
    elif interaction == 'photopion':
        s = []
        cross_section = []

        with open(f"{XSECTIONS_DIR}/xsecs_photopion_proton_sophia.txt", 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                
                values = line.strip().split(',')
                s.append(float(values[0]))  
                cross_section.append(float(values[1]))

        s = np.array(s) * GeV_to_eV**2
        cross_section = A * np.array(cross_section) * mbarn_to_m2
        eps = (s - mp**2) / (2*mp)

        return eps, cross_section
    
    else:
        raise ValueError(f"Invalid interaction type: {interaction}. Expected 'photodisintegration' or 'photopion'.")

# ----------------------------------------------------------------------------------------------------
def compute_timescales(A, Z, Gmms, interaction):

    timescales = []

    eps, cross_section = get_eps_and_cross_sections(A, Z, interaction)

    for Gmm in Gmms:
        integrand_interaction_rate = c / (2 * Gmm**2) * eps * cross_section * inelasticity(eps) * I(eps, Gmm)
        interaction_rate = simps(integrand_interaction_rate, eps)
        timescales.append(interaction_rate**-1 * s_to_yr)

    return np.array(timescales)

# ----------------------------------------------------------------------------------------------------
def write_timescales(A, Z, interaction):

    Gmms = np.logspace(8, 12, num = 100) / A
    E = Gmms * A * mp
    timescales = compute_timescales(A, Z, Gmms, interaction)   
    np.savetxt(f"{RESULTS_DIR}/timescales_{interaction}_{PARTICLES[iZ(Z)]}.dat", np.column_stack((E, timescales)), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_timescales(1, 1, 'photopion')
    # write_timescales(14, 7, 'photopion')
    # write_timescales(28, 14, 'photopion')
    # write_timescales(56, 26, 'photopion')

    # write_timescales(14, 7, 'photodisintegration')
    # write_timescales(28, 14, 'photodisintegration')
    # write_timescales(56, 26, 'photodisintegration')

# ----------------------------------------------------------------------------------------------------