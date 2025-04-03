from scipy.integrate import trapz, simps
from scipy.interpolate import interp1d
import numpy as np

RESULTS_DIR = "../results"

ls = ['gmm', 'pos', 'anu_mu', 'nu_mu', 'nu_e', 'e', 'anu_e']

c =  299792458e2    # cm/s
h = 4.135667696e-15 # eV.s
kB = 8.6173303e-5   # eV/K
mp = 0.9383e9       # eV 

eta0 = 0.313

T_IR = 3.5e-3 / kB    # K
T_OPT = 332.5e-3 / kB # K
NORM_IR = 0.15268288372409347
NORM_OPT = 2.810860522825165e-09

# ----------------------------------------------------------------------------------------------------
def photon_density(eps): # 1 / cm^-3 / eV 

    return 8 * np.pi / (h * c)**3 * eps**2 * (NORM_IR / (np.exp(eps / (kB * T_IR)) - 1) + NORM_OPT / (np.exp(eps / (kB * T_OPT)) - 1))

# ----------------------------------------------------------------------------------------------------
def get_eta_xs_and_phi_l(l):

    if l == 'gmm':
        eta = np.loadtxt("KA08_TableI.dat")[:,0] * eta0

    elif l == 'pos' or l == 'anu_mu' or l == 'nu_mu' or l == 'nu_e':
        eta = np.loadtxt("KA08_TableII.dat")[:,0] * eta0

    elif l == 'e' or l == 'anu_e':
        eta = np.loadtxt("KA08_TableIII.dat")[:,0] * eta0

    xs = np.loadtxt(f"{RESULTS_DIR}/KA08_phi_{l}.dat")[:,0]
    phi_l = np.loadtxt(f"{RESULTS_DIR}/KA08_phi_{l}.dat")[:,1:] # cm^3 / s

    return eta, xs, phi_l

# ----------------------------------------------------------------------------------------------------
def compute_mono_spectrum(l, E): # Equation (42) for a given energy E

    eta, xs, phi_l = get_eta_xs_and_phi_l(l)

    spec_mono = [] # 1 / s

    for ix in range(len(xs)):

        Ep = E / xs[ix]

        f = interp1d(eta, phi_l[ix,:], kind = 'linear', bounds_error = False, fill_value = 'extrapolate')
        eta_new = np.logspace(np.log10(eta[0]), np.log10(eta[-1]), num = 100)
        phi_l_new = f(eta_new)

        eps = eta_new * mp**2 / (4 * Ep)
        integrand = photon_density(eps) * phi_l_new

        if simps(integrand, eps) < 0:
            spec_mono.append(trapz(integrand, eps))
        else:
            spec_mono.append(simps(integrand, eps)) 

    return np.array(spec_mono)

# ----------------------------------------------------------------------------------------------------
def compute_spectrum(l, E):
        
    _, xs, _ = get_eta_xs_and_phi_l(l)

    data = np.loadtxt(f"{RESULTS_DIR}/transport_sol_1H.dat")
    Ep = data[:,0]
    fp = data[:,1]
    fp_interp = interp1d(Ep, fp, bounds_error = False, fill_value = 0)

    return simps(fp_interp(E / xs) * compute_mono_spectrum(l, E) / xs, xs)
     
# ----------------------------------------------------------------------------------------------------
def write_spectrum(l):

    Es = np.logspace(17, 21, num = 100)
    spec = []

    for E in Es:
        spec.append(compute_spectrum(l, E))

    np.savetxt(f"{RESULTS_DIR}/KA08_spectrum_{l}.dat", np.column_stack((Es, np.array(spec))), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for l in ls:
        write_spectrum(l)

# ----------------------------------------------------------------------------------------------------
