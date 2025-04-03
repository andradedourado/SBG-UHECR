from scipy.integrate import simps
from scipy.interpolate import interp1d
import numpy as np

RESULTS_DIR = "../results"

ls = ['gmm', 'nu_mu', 'e']

eV_to_TeV = 1e-12
TeV_to_eV = 1e12
mb_to_cm2 = 1e-27 

c = 2.998e10 # cm/s
nH = 125 # cm^-3

# ----------------------------------------------------------------------------------------------------
def xsection_inel(Ep): 

    Ep = Ep * eV_to_TeV
    L = np.log(Ep) # Ep in TeV

    return (34.3 + 1.88*L + 0.25*L**2) * mb_to_cm2 # cm^2

# ----------------------------------------------------------------------------------------------------
def F_l(x, Ep, l):
    
    Ep = Ep * eV_to_TeV
    L = np.log(Ep) # Ep in TeV 

    if l == 'gmm':
        
        B_gmm = 1.30 + 0.14*L + 0.011*L**2
        beta_gmm = 1 / (1.79 + 0.11*L + 0.008*L**2)
        k_gmm = 1 / (0.801 + 0.049*L + 0.014*L**2)

        F_gmm = B_gmm * np.log(x) / x \
            * ((1 - x**beta_gmm) / (1 + k_gmm * x**beta_gmm * (1 - x**beta_gmm)))**4 \
            * (1 / np.log(x) - 4 * beta_gmm * x**beta_gmm / (1 - x**beta_gmm) - 4 * k_gmm * beta_gmm * x*beta_gmm * (1 - 2 * x**beta_gmm) / (1 + k_gmm * x**beta_gmm * (1 - x**beta_gmm)))  

        return F_gmm
    
    if l == 'e':

        B_e = 1 / (69.5 + 2.65*L + 0.3*L**2)
        beta_e = 1 / (0.201 + 0.062*L + 0.00042*L**2)**(1/4)
        k_e = (0.279 + 0.141*L + 0.0172*L**2) / (0.3 + (2.3 + L)**2)

        F_e = B_e * (1 + k_e * np.log(x)**2)**3 / (x * (1 + 0.3/x**beta_e)) * (-np.log(x))**5

        return F_e
    
    if l == 'nu_mu':

        B_nu_mu_2 = 1 / (69.5 + 2.65*L + 0.3*L**2)
        beta_nu_mu_2 = 1 / (0.201 + 0.062*L + 0.00042*L**2)**(1/4)
        k_nu_mu_2 = (0.279 + 0.141*L + 0.0172*L**2) / (0.3 + (2.3 + L)**2)

        F_nu_mu_2 = B_nu_mu_2 * (1 + k_nu_mu_2 * np.log(x)**2)**3 / (x * (1 + 0.3/x**beta_nu_mu_2)) * (-np.log(x))**5

        mask = x < 0.427
        y = np.zeros_like(x)
        F_nu_mu_1 = np.zeros_like(x)

        y[mask] = x[mask] / 0.427

        B_prime = 1.75 + 0.204*L + 0.010*L**2
        beta_prime = 1 / (1.67 + 0.111*L + 0.0038*L**2)
        k_prime = 1.07 - 0.086*L + 0.002*L**2

        F_nu_mu_1[mask] = B_prime * np.log(y[mask]) / y[mask] \
                * ((1 - y[mask]**beta_prime) / (1 + k_prime * y[mask]**beta_prime * (1 - y[mask]**beta_prime)))**4 \
                * (1 / np.log(y[mask]) - 4 * beta_prime * y[mask]**beta_prime / (1 - y[mask]**beta_prime) - 4 * k_prime * beta_prime * y[mask]**beta_prime * (1 - 2 * y[mask]**beta_prime) / (1 + k_prime * y[mask]**beta_prime * (1 - y[mask]**beta_prime)))

        return F_nu_mu_1 + F_nu_mu_2

# ----------------------------------------------------------------------------------------------------
def spectrum(E, l):

    data = np.loadtxt(f"{RESULTS_DIR}/transport_sol_1H.dat")
    Ep = data[:,0]
    Jp = data[:,1]

    x = E / Ep
    integrand = xsection_inel(E/x) * Jp * F_l(x, E/x, l) / x
    
    x_inv = x[::-1]
    integrand_inv = integrand[::-1]

    mask = (x_inv >= 1e-3) & (x_inv <= 1)

    if not np.any(mask):
        return 0
    
    return c * nH * simps(integrand_inv[mask], x_inv[mask])

# ----------------------------------------------------------------------------------------------------
def write_spectrum(l):

    Es = np.logspace(-1, 3, num = 100) * TeV_to_eV  
    spec = []

    for E in Es:
        spec.append(spectrum(E, l))

    np.savetxt(f"{RESULTS_DIR}/KAB06_spectrum_{l}.dat", np.column_stack((Es, spec)), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for l in ls:
        write_spectrum(l)

# ----------------------------------------------------------------------------------------------------