from scipy.interpolate import interp1d
import numpy as np

RESULTS_DIR = "../results"
TIMESCALES_DIR = "../../timescales/results"

PARTICLES = ['1H', '4He', '14N', '28Si', '56Fe']
ZS = [1, 2, 7, 14, 26]

yr_to_s = 60 * 60 * 24 * 365.25

# ----------------------------------------------------------------------------------------------------
def iZ(Z):

    try:
        return ZS.index(Z)
    except ValueError:
        raise ValueError(f"Z ({Z}) not found in ZS.")

# ----------------------------------------------------------------------------------------------------
def injection_term(E, Z):

    return E**-2 # Implement the luminosity 

# ----------------------------------------------------------------------------------------------------
def compute_effective_timescale_from_interpolation(E, files, energy_conversion):

    interpolated_tau = []

    for file in files:
        data = np.loadtxt(f"{TIMESCALES_DIR}/{file}")
        E_table = data[:,0] * energy_conversion[file]
        tau_table = data[:,1] * yr_to_s
        tau_table = np.where(np.isinf(tau_table), np.finfo(float).max, tau_table)

        f_interp = interp1d(E_table, tau_table, kind = 'linear', bounds_error = False, fill_value = 'extrapolate')
        interpolated_tau.append(f_interp(E))

    total_tau_inv = sum(1 / tau for tau in interpolated_tau)

    return np.where(total_tau_inv > 0, 1 / total_tau_inv, np.inf)

# ----------------------------------------------------------------------------------------------------
def escaping_timescales(E, Z):

    files = ["timescales_advection.dat", f"timescales_diff_{PARTICLES[iZ(Z)]}.dat"]
    energy_conversion = {"timescales_advection.dat": 1.e18, f"timescales_diff_{PARTICLES[iZ(Z)]}.dat": 1.e18}

    return compute_effective_timescale_from_interpolation(E, files, energy_conversion)
    
# ----------------------------------------------------------------------------------------------------
def energy_loss_timescales(E, Z):

    files = [f"timescales_spal_{PARTICLES[iZ(Z)]}.dat", f"timescales_pairproduction_{PARTICLES[iZ(Z)]}.dat",
            f"timescales_photopion_{PARTICLES[iZ(Z)]}.dat", f"timescales_photodisintegration_{PARTICLES[iZ(Z)]}.dat"]
    
    energy_conversion = {f"timescales_spal_{PARTICLES[iZ(Z)]}.dat": 1.e18, f"timescales_pairproduction_{PARTICLES[iZ(Z)]}.dat": 1.0, 
                        f"timescales_photopion_{PARTICLES[iZ(Z)]}.dat": 1.0, f"timescales_photodisintegration_{PARTICLES[iZ(Z)]}.dat": 1.0}

    return compute_effective_timescale_from_interpolation(E, files, energy_conversion)

# ----------------------------------------------------------------------------------------------------
def total_timescales(E, Z):

    return (1 / escaping_timescales(E, Z) + 1 / energy_loss_timescales(E, Z))**-1

# ----------------------------------------------------------------------------------------------------
def cr_equilibrium_density(E, Z):

    return injection_term(E, Z) * total_timescales(E, Z)

# ----------------------------------------------------------------------------------------------------
def cr_escaping_emissivity(E, Z):

    return cr_equilibrium_density(E, Z) / escaping_timescales(E, Z)

# ----------------------------------------------------------------------------------------------------
def write_cr_escaping_emissivity(E, Z):

    np.savetxt(f"{RESULTS_DIR}/cr_escaping_emissivity_{PARTICLES[iZ(Z)]}.dat", np.column_stack((E, cr_escaping_emissivity(E, Z))), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    E = np.logspace(16, 21, num = 100)

    write_cr_escaping_emissivity(E, 26)

# ----------------------------------------------------------------------------------------------------