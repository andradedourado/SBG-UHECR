from scipy.integrate import quad
from scipy.integrate import trapz, simps
from scipy.interpolate import interp1d
import numpy as np

RESULTS_DIR = "../results"

Ecuts = np.array([0.1, 1, 10, 1000]) * 3e20 # eV
ls = ['gmm', 'pos', 'anu_mu', 'nu_mu', 'nu_e', 'e', 'anu_e']

erg_to_eV = 6.242e11

c =  299792458e2    # cm/s
h = 4.135667696e-15 # eV.s
kB = 8.6173303e-5   # eV/K
mp = 0.9383e9       # eV 

eta0 = 0.313

# ----------------------------------------------------------------------------------------------------
def fp_norm_integrand(Ep, Ecut): # cm^-3

    return Ep * Ep**-2 * np.exp(-Ep/Ecut)

# ----------------------------------------------------------------------------------------------------
def fp(Ep, Ecut):

    if Ecut == 0.1 * 3e20:
        A = erg_to_eV / quad(fp_norm_integrand, 1e9, 1e23, args = (Ecut))[0] 
    else:
        A = erg_to_eV / quad(fp_norm_integrand, 1e9, 1e24, args = (Ecut))[0] 
    
    return A * Ep**-2 * np.exp(-Ep/Ecut)

# ----------------------------------------------------------------------------------------------------
def photon_density_CMB(eps, T = 2.7): # 1 / cm^-3 / eV 

    return 8 * np.pi * eps**2 / (h * c)**3 / (np.exp(eps / (kB * T)) - 1)

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
def write_mono_spectrum(l, Ep): # Equation (42)

    eta, xs, phi_l = get_eta_xs_and_phi_l(l)

    spec = [] # 1 / s

    for ix in range(len(xs)):

        f = interp1d(eta, phi_l[ix,:], kind = 'linear', bounds_error = False, fill_value = 'extrapolate')
        eta_new = np.logspace(np.log10(eta[0]), np.log10(eta[-1]), num = 100)
        phi_l_new = f(eta_new)

        eps = eta_new * mp**2 / (4 * Ep)
        integrand = photon_density_CMB(eps) * phi_l_new

        if simps(integrand, eps) < 0:
            spec.append(trapz(integrand, eps))
        else:
            spec.append(simps(integrand, eps)) 

    np.savetxt(f"{RESULTS_DIR}/KA08_mono_spectrum_CMB_Ep{int(np.log10(Ep)):02d}_{l}.dat", np.column_stack((xs, np.array(spec))), fmt = "%.15e")

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
        integrand = photon_density_CMB(eps) * phi_l_new

        if simps(integrand, eps) < 0:
            spec_mono.append(trapz(integrand, eps))
        else:
            spec_mono.append(simps(integrand, eps)) 

    return np.array(spec_mono)

# ----------------------------------------------------------------------------------------------------
def write_spectrum(l, Ecut):

    _, xs, _ = get_eta_xs_and_phi_l(l)

    Es = np.logspace(17, 21, num = 100)
    
    spec = []

    for E in Es:
        integrand = fp(E / xs, Ecut) * compute_mono_spectrum(l, E) / xs
        spec.append(simps(integrand, xs))

    Ecut_str = f"{int(Ecut / 3e20)}" if (Ecut / 3e20).is_integer() else f"{Ecut / 3e20:.1f}".replace(".", "_")
    np.savetxt(f"{RESULTS_DIR}/KA08_spectrum_CMB_{Ecut_str}_{l}.dat", np.column_stack((Es, np.array(spec))), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for l in ls:
        # write_mono_spectrum(l, 1e20)
        # write_mono_spectrum(l, 1e21)
        for Ecut in Ecuts:
            write_spectrum(l, Ecut)

# ----------------------------------------------------------------------------------------------------
