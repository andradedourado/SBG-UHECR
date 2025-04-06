import numpy as np

RESULTS_DIR = "../results"

INTERACTIONS = ['advection', 'diff', 'MR19_diff', 'spal']
PARTICLES = ['1H', '4He', '14N', '28Si', '56Fe']
ZSS = [1, 2, 7, 14, 26]

EeV_to_TeV = 1.e6
TeV_to_EeV = 1.e-6
km_to_cm = 100000   
pc_to_cm = 3.0857e18 
mb_to_cm2 = 1e-27  
s_to_yr = 1 / (60 * 60 * 24 * 365.25)

R = 225 * pc_to_cm # cm 
B = 200 # μG
vW = 500 * km_to_cm # cm/s
nISM = 125 # cm^-3

c = 3e10 # cm/s

# ----------------------------------------------------------------------------------------------------
def iZs(Zs):

    try:
        return ZSS.index(Zs)
    except ValueError:
        raise ValueError(f"Zs ({Zs}) not found in ZSS.")
    
# ----------------------------------------------------------------------------------------------------
def diffusion_coefficient(E, Z, dlt = 5/3, lc = 1):

    rL = 1.081e3 / Z * (E / B) # pc 

    D = np.zeros_like(E)

    mask_low = rL < lc
    mask_high = ~mask_low  # rL >= lc

    D[mask_low] = c * (rL[mask_low] * pc_to_cm) ** (2 - dlt) * (lc * pc_to_cm) ** (dlt - 1) / 3
    D[mask_high] = c * (lc * pc_to_cm) / 3 * (rL[mask_high] / lc) ** 2

    return D

# ----------------------------------------------------------------------------------------------------
def MR19_diffusion_coefficient(E, Z, m = 5/3, lc = 1):

    aI = 0.9
    aL = 0.23
    Ec = 0.9e-3 * Z * B * lc # EeV

    return c / 3 * (lc * pc_to_cm) * (4 * (E / Ec)**2 + aI * (E / Ec) + aL * (E / Ec)**(2 - m)) # cm^2 / s

# ----------------------------------------------------------------------------------------------------
def spallation_cross_section(E, A):

    Eth = 1.22e-3 * TeV_to_EeV # EeV 
    L = np.log(E * EeV_to_TeV)
    
    return A**(2/3) * (34.3 + 1.88 * L + 0.25 * L**2) * (1 - (Eth/E)**4)**2 * mb_to_cm2 # cm^2

# ----------------------------------------------------------------------------------------------------
def advection_time(R, vW):

    return R / vW * s_to_yr

# ----------------------------------------------------------------------------------------------------
def diffusion_times(E, Z):

    diffusion_times = R**2 / diffusion_coefficient(E, Z)
    min_time = R / c
    return np.where(diffusion_times > min_time, diffusion_times, min_time) * s_to_yr

# ----------------------------------------------------------------------------------------------------
def MH19_diffusion_times(E, Z):

    diffusion_times = R**2 / MR19_diffusion_coefficient(E, Z)
    min_time = R / c
    return np.where(diffusion_times > min_time, diffusion_times, min_time) * s_to_yr

# ----------------------------------------------------------------------------------------------------
def spallation_times(E, A):

    return (0.5 * nISM * spallation_cross_section(E, A) * c)**-1 * s_to_yr

# ----------------------------------------------------------------------------------------------------
def compute_timescales(interaction, E, Z, A):

    if interaction == 'advection':
        return advection_time(R, vW)

    elif interaction == 'diff':
        return diffusion_times(E, Z)

    elif interaction == 'MR19_diff':
        return MH19_diffusion_times(E, Z)

    elif interaction == 'spal':
        return spallation_times(E, A)

    else:
        raise ValueError(f"Error: The interaction type '{interaction}' is not recognized. Please choose from 'advection', 'diff', 'photo', or 'spal'.")

# ----------------------------------------------------------------------------------------------------
def write_timescales(interaction, Z, A, nbins = 100):

    E_values = np.logspace(-2, 3, num = nbins)
    timescales = compute_timescales(interaction, E_values, Z, A)

    if interaction == 'advection':
        with open(f"{RESULTS_DIR}/timescales_{interaction}.dat", 'w') as f:
            for E in E_values:
                f.write(f'{E:.15e}\t{compute_timescales(interaction, E_values, Z, A):.15e}\n')

    else:
        with open(f"{RESULTS_DIR}/timescales_{interaction}_{PARTICLES[iZs(Z)]}.dat", 'w') as f:
            for E, timescale in zip(E_values, timescales):
                f.write(f'{E:.15e}\t{timescale:.15e}\n')

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_timescales('advection', 0, 0)
    write_timescales('diff', 1, 1)
    write_timescales('diff', 26, 56)
    write_timescales('MR19_diff', 1, 1)
    write_timescales('MR19_diff', 26, 56)
    write_timescales('spal', 1, 1)
    write_timescales('spal', 26, 56)

# ----------------------------------------------------------------------------------------------------