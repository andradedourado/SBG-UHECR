import numpy as np 

RESULTS_DIR = "../results"

ls = ['gmm', 'pos', 'anu_mu', 'nu_mu', 'nu_e', 'e', 'anu_e']

r = 0.146 # m_pi / m_p 

# ----------------------------------------------------------------------------------------------------
def heaviside(x, x0):

    return (np.array(x) >= x0).astype(int)

# ----------------------------------------------------------------------------------------------------
def phi_gmm(x):

    KA08_TableI = np.loadtxt(f"KA08_TableI.dat")

    eta_over_eta0 = KA08_TableI[:,0]
    s_gmm = KA08_TableI[:,1]
    dlt_gmm = KA08_TableI[:,2]
    B_gmm = KA08_TableI[:,3] * 10**KA08_TableI[:,4] # cm^3 / s

    eta = eta_over_eta0 * 0.313 # Equation (16)
    x_minus = 1 / (2 * (1 + eta)) * (eta + r**2 - np.sqrt((eta - r**2 - 2*r) * (eta - r**2 + 2*r)))  
    x_plus = 1 / (2 * (1 + eta)) * (eta + r**2 + np.sqrt((eta - r**2 - 2*r) * (eta - r**2 + 2*r))) 

    phi = np.zeros_like(eta_over_eta0)
    y = np.zeros_like(eta_over_eta0)

    mask_x_between_xlimits = (x_minus < x) & (x_plus > x)    
    mask_x_lt_xminus = x < x_minus
    mask_x_ge_xplus = x >= x_plus

    y[mask_x_between_xlimits] = (x - x_minus[mask_x_between_xlimits]) / (x_plus[mask_x_between_xlimits] - x_minus[mask_x_between_xlimits])
    phi[mask_x_between_xlimits] = B_gmm[mask_x_between_xlimits] * np.exp(-s_gmm[mask_x_between_xlimits] * (np.log(x / x_minus[mask_x_between_xlimits]))**dlt_gmm[mask_x_between_xlimits]) * (np.log(2 / (1 + y[mask_x_between_xlimits]**2)))**(2.5 + 0.4 * np.log(eta_over_eta0[mask_x_between_xlimits]))
    phi[mask_x_lt_xminus] = B_gmm[mask_x_lt_xminus] * np.log(2)**(2.5 + 0.4 * np.log(eta_over_eta0[mask_x_lt_xminus]))
    phi[mask_x_ge_xplus] = 0

    return phi

# ----------------------------------------------------------------------------------------------------
def phi_l(x, l):

    KA08_TableII = np.loadtxt(f"KA08_TableII.dat")
    KA08_TableIII = np.loadtxt(f"KA08_TableIII.dat") 

    if l == 'pos':

        eta_over_eta0 = KA08_TableII[:,0]
        s_l = KA08_TableII[:,1]
        dlt_l = KA08_TableII[:,2]
        B_l = KA08_TableII[:,3] * 10**KA08_TableII[:,4] # cm^3 / s 

    elif l == 'anu_mu':

        eta_over_eta0 = KA08_TableII[:,0]
        s_l = KA08_TableII[:,5]
        dlt_l = KA08_TableII[:,6]
        B_l = KA08_TableII[:,7] * 10**KA08_TableII[:,8] # cm^3 / s

    elif l == 'nu_mu':

        eta_over_eta0 = KA08_TableII[:,0]
        s_l = KA08_TableII[:,9]
        dlt_l = KA08_TableII[:,10]
        B_l = KA08_TableII[:,11] * 10**KA08_TableII[:,12] # cm^3 / s

    elif l == 'nu_e':

        eta_over_eta0 = KA08_TableII[:,0]
        s_l = KA08_TableII[:,13]
        dlt_l = KA08_TableII[:,14]
        B_l = KA08_TableII[:,15] * 10**KA08_TableII[:,16] # cm^3 / s

    elif l == 'e':

        eta_over_eta0 = KA08_TableIII[:,0]
        s_l = KA08_TableIII[:,1]
        dlt_l = KA08_TableIII[:,2]
        B_l = KA08_TableIII[:,3] * 10**KA08_TableIII[:,4] # cm^3 / s

    elif l == 'anu_e':
        
        eta_over_eta0 = KA08_TableIII[:,0]
        s_l = KA08_TableIII[:,5]
        dlt_l = KA08_TableIII[:,6]
        B_l = KA08_TableIII[:,7] * 10**KA08_TableIII[:,8] # cm^3 / s

    eta = eta_over_eta0 * 0.313 # Equation (16)
    x_minus = 1 / (2 * (1 + eta)) * (eta + r**2 - np.sqrt((eta - r**2 - 2*r) * (eta - r**2 + 2*r)))  
    x_plus = 1 / (2 * (1 + eta)) * (eta + r**2 + np.sqrt((eta - r**2 - 2*r) * (eta - r**2 + 2*r))) 

    if l == 'pos' or l == 'anu_mu' or l == 'nu_e':
        
        x_prime_minus = x_minus / 4
        x_prime_plus = x_plus
        psi = 2.5 + 1.4 * np.log(eta_over_eta0)

    elif l == 'nu_mu':
        
        x_prime_minus = 0.427 * x_minus
        x_prime_plus = np.zeros_like(x_plus)

        mask_rho_lt_2_14 = eta_over_eta0 < 2.14
        mask_rho_between_2_14_10 = (eta_over_eta0 > 2.14) & (eta_over_eta0 < 10)
        mask_rho_gt_10 = eta_over_eta0 > 10

        x_prime_plus[mask_rho_lt_2_14] = 0.427 * x_plus[mask_rho_lt_2_14]
        x_prime_plus[mask_rho_between_2_14_10] = (0.427 + 0.0729 * (eta_over_eta0[mask_rho_between_2_14_10] - 2.14)) * x_plus[mask_rho_between_2_14_10]
        x_prime_plus[mask_rho_gt_10] = x_plus[mask_rho_gt_10]
        
        psi = 2.5 + 1.4 * np.log(eta_over_eta0)

    elif l == 'e' or l == 'anu_e':

        x_min = 1 / (2 * (1 + eta)) * (eta - 2*r - np.sqrt(eta * (eta - 4 * r * (1 + r))))
        x_max = 1 / (2 * (1 + eta)) * (eta - 2*r + np.sqrt(eta * (eta - 4 * r * (1 + r))))
        
        x_prime_minus = x_min / 2
        x_prime_plus = x_max
        
        psi = 6 * (1 - np.exp(1.5 * (4 - eta_over_eta0))) * heaviside(eta_over_eta0, 4)

    phi = np.zeros_like(eta_over_eta0)
    y_prime = np.zeros_like(eta_over_eta0)

    mask_x_between_xprime = (x > x_prime_minus) & (x < x_prime_plus)
    mask_x_lt_xprime_minus = x < x_prime_minus
    mask_x_ge_xprime_plus = x >= x_prime_plus

    y_prime[mask_x_between_xprime] = (x - x_prime_minus[mask_x_between_xprime]) / (x_prime_plus[mask_x_between_xprime] - x_prime_minus[mask_x_between_xprime])
    phi[mask_x_between_xprime] = B_l[mask_x_between_xprime] * np.exp(-s_l[mask_x_between_xprime] * (np.log(x / x_prime_minus[mask_x_between_xprime]))**dlt_l[mask_x_between_xprime]) * (np.log(2 / (1 + y_prime[mask_x_between_xprime]**2)))**psi[mask_x_between_xprime]
    phi[mask_x_lt_xprime_minus] = B_l[mask_x_lt_xprime_minus] * (np.log(2))**psi[mask_x_lt_xprime_minus]
    phi[mask_x_ge_xprime_plus] = 0 

    return phi

# ----------------------------------------------------------------------------------------------------
def write_phi(l):

    phi = []
    xs = np.logspace(-4, 1, num = 1000)

    for x in xs:
        if l == 'gmm':
            phi.append(phi_gmm(x))
        else:
            phi.append(phi_l(x, l))
    
    phi = np.array(phi)
    np.savetxt(f"{RESULTS_DIR}/KA08_phi_{l}.dat", np.column_stack((xs, phi)), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for l in ls:
        write_phi(l)

# ----------------------------------------------------------------------------------------------------