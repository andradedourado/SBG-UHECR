from scipy.integrate import quad
import numpy as np

RESULTS_DIR = "../results"

alps = [2, 1.5]
ls = ['gmm', 'nu_mu', 'e']

erg_to_TeV = 0.6242
mb_to_cm2 = 1e-27 

c = 2.998e10 # cm/s
nH = 1 # cm^-3

# ----------------------------------------------------------------------------------------------------
def xsection_inel(Ep): 

    L = np.log(Ep) # Ep in TeV

    return (34.3 + 1.88*L + 0.25*L**2) * mb_to_cm2 # cm^2

# ----------------------------------------------------------------------------------------------------
def Jp_norm_integrand(Ep, alp, beta = 1, E0 = 1000): # cm^-3

    return Ep / Ep**alp * np.exp(-(Ep/E0)**beta)

# ----------------------------------------------------------------------------------------------------
def Jp(Ep, alp, beta = 1, E0 = 1000): # cm^-3 TeV^-1

    A = erg_to_TeV / quad(Jp_norm_integrand, 1, np.inf, args = (alp))[0] 
    return A / Ep**alp * np.exp(-(Ep/E0)**beta) 

# ----------------------------------------------------------------------------------------------------
def F_l(x, Ep, l):

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

        if x < 0.427:

            y = x / 0.427
            B_prime = 1.75 + 0.204*L + 0.010*L**2
            beta_prime = 1 / (1.67 + 0.111*L + 0.0038*L**2)
            k_prime = 1.07 - 0.086*L + 0.002*L**2

            F_nu_mu_1 = B_prime * np.log(y) / y \
                * ((1 - y**beta_prime) / (1 + k_prime * y**beta_prime * (1 - y**beta_prime)))**4 \
                * (1 / np.log(y) - 4 * beta_prime * y**beta_prime / (1- y**beta_prime) - 4 * k_prime * beta_prime * y**beta_prime * (1 - 2 * y**beta_prime) / (1 + k_prime * y**beta_prime * (1 - y**beta_prime))) 

        elif x >= 0.427:

            F_nu_mu_1 = 0

        return F_nu_mu_1 + F_nu_mu_2

# ----------------------------------------------------------------------------------------------------
def spectrum_integrand(x, E, alp, l):

    return xsection_inel(E/x) * Jp(E/x, alp) * F_l(x, E/x, l) / x

# ----------------------------------------------------------------------------------------------------
def spectrum(E, alp, l):

    return c * nH * quad(spectrum_integrand, 1e-3, 1, args = (E, alp, l))[0]

# ----------------------------------------------------------------------------------------------------
def write_spectrum(l, alp):

    Es = np.logspace(-1, 3, num = 100) # TeV 
    spec = []

    for E in Es:
        spec.append(spectrum(E, alp, l))

    if alp == 2:
        np.savetxt(f"{RESULTS_DIR}/KAB06_spectrum_2_{l}.dat", np.column_stack((Es, spec)), fmt = "%.15e")
    
    elif alp == 1.5:
        np.savetxt(f"{RESULTS_DIR}/KAB06_spectrum_1_5_{l}.dat", np.column_stack((Es, spec)), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for l in ls:
        for alp in alps:
            write_spectrum(l, alp)

# ----------------------------------------------------------------------------------------------------


