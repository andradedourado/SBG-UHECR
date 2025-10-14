from scipy.integrate import quad
import numpy as np 

RESULTS_DIR = "../results"
SIMULATIONS_RESULTS_DIR = "../../simulations/results/"

ENERGY_EDGES = np.logspace(0, 4, num = 80)
ENERGY_BIN_CENETRS = np.logspace(4./(2.*(len(ENERGY_EDGES) - 1)), 4. - 4./(2.*(len(ENERGY_EDGES) - 1)), num = len(ENERGY_EDGES) - 1)
PARTICLES = ['1H', '4He', '14N', '28Si', '56Fe']
ZSS = [1, 2, 7, 14, 26]

galaxy_data = np.genfromtxt("../../galaxies/scripts/starburst_galaxies.dat", dtype = None, encoding = None)

eV_to_erg = 1.60218e-12

Emin = 10**17.8 # eV

E0 = 1e18

# Implement the emissivity from the leaky-box model
# Learn how to normalize the outputs of the simulations

# ----------------------------------------------------------------------------------------------------
def iZs(Zs):

    try:
        return ZSS.index(Zs)
    except ValueError:
        raise ValueError(f"Zs ({Zs}) not found in ZSS.")

# ----------------------------------------------------------------------------------------------------
def w_sim(Es):

    return Es * 1e18 / E0

# ----------------------------------------------------------------------------------------------------
def w_spec(Es, Zs):

    L0, Gmm, Rcut = generation_rate_parameters(Zs)

    Es = Es * 1e18

    mask_low = Es <= Zs * Rcut
    mask_high = ~mask_low

    w_spec = np.zeros_like(Es)

    w_spec[mask_low] = (Es[mask_low] / E0)**-Gmm
    w_spec[mask_high] = (Es[mask_high] / E0)**-Gmm * np.exp(1 - Es[mask_high] / (Zs * Rcut))

    if Zs == 1:
        return w_spec * L0 / (quad(integrand_w_spec, Emin, 1e23, args = (Zs))[0] * eV_to_erg**2)
    else:
        return w_spec * [0.0, 0.245, 0.681, 0.049, 0.025][iZs(Zs)] * L0 / (quad(integrand_w_spec, Emin, 1e23, args = (Zs))[0] * eV_to_erg**2)

# ----------------------------------------------------------------------------------------------------
def integrand_w_spec(Es, Zs):

    _, Gmm, Rcut = generation_rate_parameters(Zs)

    if Es <= Zs * Rcut: 
        return Es * (Es / E0)**-Gmm
    elif Es > Zs * Rcut: 
        return Es * (Es / E0)**-Gmm * np.exp(1 - Es / (Zs * Rcut))
    
# ----------------------------------------------------------------------------------------------------
def generation_rate_parameters(Z):

    if Z == 1:
        L0 = 6.54e44 # erg Mpc^-3 yr^-1
        Gmm = 3.34
        Rcut = 10**19.3 # V
        return L0, Gmm, Rcut

    else:
        L0 = 5e44 # erg Mpc^-3 yr^-1
        Gmm = -1.47
        Rcut = 10**18.19 # V
        return L0, Gmm, Rcut

# ----------------------------------------------------------------------------------------------------
def compute_intensity(Zs):

    spec = np.zeros_like(ENERGY_BIN_CENETRS)

    for igal in range(len(galaxy_data)):
        if galaxy_data[igal][2] == 'NGC4038/9':
            data = np.loadtxt(f"{SIMULATIONS_RESULTS_DIR}/{PARTICLES[iZs(Zs)]}/events_NGC4038_9.txt")
        else:
            data = np.loadtxt(f"{SIMULATIONS_RESULTS_DIR}/{PARTICLES[iZs(Zs)]}/events_{galaxy_data[igal][2]}.txt")
    
        spec += np.histogram(data[:,2], bins = ENERGY_EDGES, weights = w_sim(data[:,4]) * w_spec(data[:,4], Zs))[0] / (4 * np.pi * galaxy_data[igal][8])**2

    return spec 

# ----------------------------------------------------------------------------------------------------
def write_intensity(Zs):

    spec = compute_intensity(Zs)
    np.savetxt(f"{RESULTS_DIR}/spec_{PARTICLES[iZs(Zs)]}.dat", np.column_stack((ENERGY_BIN_CENETRS * 1e18, spec / (ENERGY_BIN_CENETRS * 1e18))), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    
    for Zs in ZSS:
        write_intensity(Zs)

# ----------------------------------------------------------------------------------------------------