from scipy.integrate import simps
import matplotlib.pyplot as plt
import numpy as np 

REFERENCES_DIR = "../references"
RESULTS_DIR = "../results"
XSECTIONS_DIR = "../../cross-sections/scripts"

cm2_to_m2 = 1e-4 
mbarn_to_m2 = 1.e-31
eV_to_GeV = 1e-9
GeV_to_eV = 1.e9
pc_to_cm = 3.0857e18
s_to_yr = 1 / (60 * 60 * 24 * 365.25)

c =  299792458 # m/s
mp = 0.9383e9  # eV 

D = 3.7e6 * pc_to_cm # cm  
R = 225 * pc_to_cm   # cm 

# ----------------------------------------------------------------------------------------------------
def photon_density_integral(eps_prime_arr, Gmm):

    integral_IR = []
    integral_OPT = []

    data_IR = np.loadtxt(f"{REFERENCES_DIR}/SED_Condorelli_IR.dat")
    data_OPT = np.loadtxt(f"{REFERENCES_DIR}/SED_Condorelli_OPT.dat")

    eps_IR = data_IR[:,0] # eV
    eps_OPT = data_OPT[:,0] # eV
    
    photon_density_IR = 9/4 * (D/R)**2 / c * data_IR[:,1] / (data_IR[:,0] * eV_to_GeV)**2 / GeV_to_eV / cm2_to_m2
    photon_density_OPT = 9/4 * (D/R)**2 / c * data_OPT[:,1] / (data_OPT[:,0] * eV_to_GeV)**2 / GeV_to_eV / cm2_to_m2

    for eps_prime in eps_prime_arr:

        mask_IR = eps_IR > eps_prime / (2 * Gmm)
        mask_OPT = eps_OPT > eps_prime / (2 * Gmm)

        if np.count_nonzero(photon_density_IR[mask_IR] != 0) >= 2:
            integral_IR.append(simps(photon_density_IR[mask_IR] / eps_IR[mask_IR]**2, eps_IR[mask_IR]))
        else:
            integral_IR.append(0)  

        if np.count_nonzero(photon_density_OPT[mask_OPT] != 0) >= 2:
            integral_OPT.append(simps(photon_density_OPT[mask_OPT] / eps_OPT[mask_OPT]**2, eps_OPT[mask_OPT]))
        else:
            integral_OPT.append(0)

    integral_IR = np.array(integral_IR)
    integral_OPT = np.array(integral_OPT)

    return integral_IR + integral_OPT

# ----------------------------------------------------------------------------------------------------
def Y(x, inelasticity): 

    Y_inf = 0.47
    xb = 6e9 # eV
    dlt = 0.33
    s = 0.15

    if inelasticity == True:
        return Y_inf * (x / xb)**dlt / (1 + (x / xb)**(dlt / s))**s
    
    elif inelasticity == False:
        return 1
    
# ----------------------------------------------------------------------------------------------------
def get_eps_prime_and_cross_sections(A):

    s = []
    cross_section = []

    with open(f"{XSECTIONS_DIR}/xsecs_photopion_proton_sophia.txt", 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
                
            values = line.strip().split(',')
            s.append(float(values[0]))  
            cross_section.append(float(values[1]))

    s = np.array(s) * GeV_to_eV**2
    cross_section = A * np.array(cross_section) * mbarn_to_m2
    eps_prime = (s - mp**2) / (2*mp) # eV

    return eps_prime, cross_section

# ----------------------------------------------------------------------------------------------------
def write_timescales_xchecks(A, inelasticity):

    Gmms = np.logspace(8, 12, num = 100) / A
    E = Gmms * A * mp

    timescales = []

    eps_prime, cross_section = get_eps_prime_and_cross_sections(A)

    for Gmm in Gmms:
        integrand_interaction_rate = c / (2 * Gmm**2) * eps_prime * cross_section * Y(eps_prime, inelasticity) * photon_density_integral(eps_prime, Gmm)
        interaction_rate = simps(integrand_interaction_rate, eps_prime)
        timescales.append(interaction_rate**-1 * s_to_yr)

    if inelasticity == True: 
        np.savetxt(f"{RESULTS_DIR}/timescales_photopion_xchecks_1H_wInelasticity.dat", np.column_stack((E, timescales)), fmt = "%.15e", delimiter = "\t")

    elif inelasticity == False: 
        np.savetxt(f"{RESULTS_DIR}/timescales_photopion_xchecks_1H_woInelasticity.dat", np.column_stack((E, timescales)), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_timescales_xchecks(1, True)
    write_timescales_xchecks(1, False)

# ----------------------------------------------------------------------------------------------------
