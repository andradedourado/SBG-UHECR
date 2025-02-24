from scipy.integrate import simps
import ast
import numpy as np 
import subprocess

RESULTS_DIR = "../results"
PARTICLES = ['14N', '28Si', '56Fe']
ZS = [7, 14, 26]

mbarn_to_m2 = 1.e-31
MeV_to_eV = 1.e6     
s_to_yr = 1 / (60 * 60 * 24 * 365.25) 

c =  299792458       # m/s
hbar = 6.5821220e-16 # eV.s
kB = 8.6173303e-5    # eV/K
mp = 0.9383e9        # 1 GeV 

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
        ['python3', 'get_cross_section_TENDL-2023.py', str(A), str(Z)],
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
def compute_photodisintegration_timescales(A, Z, Gmm):
         
    eps, cross_section = execute_get_cross_section_TENDL2023(A, Z)
    mask = cross_section > 0
    eps = eps[mask] * MeV_to_eV 
    cross_section = cross_section[mask] * mbarn_to_m2

    integrand_interaction_rate = c / (2 * Gmm**2) * eps * cross_section * I(eps, Gmm)
    interaction_rate = simps(integrand_interaction_rate, eps)
    
    return interaction_rate**-1 * s_to_yr

# ----------------------------------------------------------------------------------------------------
def write_photodisintegration_timescales(A, Z):

    Gmm = np.logspace(8, 12, num = 100) / A
    E = Gmm * A * mp
    timescales = np.zeros_like(Gmm)

    f = open(f"{RESULTS_DIR}/timescales_photodisintegration_{PARTICLES[iZ(Z)]}.dat", 'w')

    for iE in range(len(E)):
        timescales[iE] = compute_photodisintegration_timescales(A, Z, Gmm[iE])
        f.write(str('{:.15e}'.format(E[iE])) + '\t')
        f.write(str('{:.15e}'.format(timescales[iE])) + '\n')
        
    f.close()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_photodisintegration_timescales(14, 7)
    write_photodisintegration_timescales(28, 14)
    write_photodisintegration_timescales(56, 26)

# ----------------------------------------------------------------------------------------------------