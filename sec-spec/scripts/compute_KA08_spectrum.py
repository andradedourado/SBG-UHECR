from scipy.integrate import simps
import numpy as np

RESULTS_DIR = "../results"

ls = ['gmm', 'pos', 'anu_mu', 'nu_mu', 'nu_e', 'e', 'anu_e']

c =  299792458e2    # cm/s
h = 4.135667696e-15 # eV.s
kB = 8.6173303e-5   # eV/K
mp = 0.9383e9       # eV 

eta0 = 0.313 

# ----------------------------------------------------------------------------------------------------
def photon_density(eps, T = 2.7): # 1 / cm^-3 / eV 

    return 8 * np.pi * eps**2 / (h * c)**3 / (np.exp(eps / (kB * T)) - 1)

# ----------------------------------------------------------------------------------------------------
def write_mono_spectrum(l, Ep): # Equation (42)

    if l == 'gmm':
        eta = np.loadtxt("KA08_TableI.dat")[:,0] * eta0

    elif l == 'pos' or l == 'anu_mu' or l == 'nu_mu' or l == 'nu_e':
        eta = np.loadtxt("KA08_TableII.dat")[:,0] * eta0

    elif l == 'e' or l == 'anu_e':
        eta = np.loadtxt("KA08_TableIII.dat")[:,0] * eta0

    xs = np.loadtxt(f"{RESULTS_DIR}/KA08_phi_{l}.dat")[:,0]
    phi_l = np.loadtxt(f"{RESULTS_DIR}/KA08_phi_{l}.dat")[:,1:] # cm^3 / s
    spec = [] # 1 / s

    eps = eta * mp**2 / (4 * Ep)

    for ix in range(len(xs)):
        integrand = photon_density(eps) * phi_l[ix,:]
        spec.append(simps(integrand, eps))

    spec = np.array(spec)
    np.savetxt(f"{RESULTS_DIR}/KA08_mono_spectrum_CMB_Ep{int(np.log10(Ep)):02d}_{l}.dat", np.column_stack((xs, spec)), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for l in ls:
        # write_mono_spectrum(l, 1e20)
        write_mono_spectrum(l, 1e21)

# ----------------------------------------------------------------------------------------------------
