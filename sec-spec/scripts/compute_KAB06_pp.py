import numpy as np 

RESULTS_DIR = "../results"

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

        y = x / 0.427
        B_prime = 1.75 * 0.204*L + 0.010*L**2
        beta_prime = 1 / (1.67 + 0.111*L + 0.0038*L**2)
        k_prime = 1.07 - 0.086*L + 0.002*L**2

        F_nu_mu_1 = B_prime * np.log(y) / y \
            * ((1 - y**beta_prime) / (1 + k_prime * y**beta_prime * (1 - y**beta_prime)))**4 \
            * (1 / np.log(y) - 4 * beta_prime * y**beta_prime / (1- y**beta_prime) - 4 * k_prime * beta_prime * y**beta_prime * (1 - 2 * y**beta_prime) / (1 + k_prime * y**beta_prime * (1 - y**beta_prime))) 

        return F_nu_mu_1 + F_nu_mu_2

# ----------------------------------------------------------------------------------------------------
def write_F(l):
    
    xs = np.logspace(-3, 0, num = 100)

    if l == 'gmm':

        F_1TeV = F_l(xs, 1, l) 
        F_30TeV = F_l(xs, 30, l) 
        F_300TeV = F_l(xs, 300, l) 
        F_3000TeV = F_l(xs, 3000, l)

        np.savetxt(f"{RESULTS_DIR}/KAB06_F_{l}.dat", np.column_stack((xs, F_1TeV, F_30TeV, F_300TeV, F_3000TeV)), fmt = "%.15e")

    elif l == 'e': 

        F_0_1TeV = F_l(xs, 0.1, l) 
        F_100TeV = F_l(xs, 100, l) 
        F_1000TeV = F_l(xs, 1000, l)

        np.savetxt(f"{RESULTS_DIR}/KAB06_F_{l}.dat", np.column_stack((xs, F_0_1TeV, F_100TeV, F_1000TeV)), fmt = "%.15e")

    elif l == 'nu_mu':
            
        mask = xs < 0.427
        xs = xs[mask]
        F_0_1TeV = F_l(xs, 0.1, l) 
        F_100TeV = F_l(xs, 100, l) 
        F_1000TeV = F_l(xs, 1000, l)

        np.savetxt(f"{RESULTS_DIR}/KAB06_F_{l}.dat", np.column_stack((xs, F_0_1TeV, F_100TeV, F_1000TeV)), fmt = "%.15e")
     
# ----------------------------------------------------------------------------------------------------
if __name__ in '__main__':

    write_F('gmm')
    write_F('e')
    write_F('nu_mu')

# ----------------------------------------------------------------------------------------------------