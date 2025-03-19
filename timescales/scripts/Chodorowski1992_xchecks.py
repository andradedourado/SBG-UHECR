from scipy.special import zeta
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'legend.fontsize': 'large',
'legend.title_fontsize': 'large',
'axes.labelsize': 'x-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

FIGURES_DIR = "../figures"
REFERENCES_DIR = "../references"

# ----------------------------------------------------------------------------------------------------
def phi_01(kappa): # kappa must be a scalar

    c1 = 0.8048
    c2 = 0.1459
    c3 = 1.137 * 1.e-3
    c4 = -3.879 * 1.e-6
    c = [c1, c2, c3, c4]

    d0 = -170 + 84*np.log(2) - 16*np.log(2)**2 + np.pi**2/3 * (10 - 4*np.log(2)) + 8*zeta(3)
    d1 = 88 - 40*np.log(2) + 8*np.log(2)**2 - 4/3*np.pi**2
    d2 = -20 + 8*np.log(2)
    d3 = 8/3
    d = [d0, d1, d2, d3]

    f1 = 2.910
    f2 = 78.35
    f3 = 1837
    f = [f1, f2, f3]

    if kappa < 25:

        sum_c_term = 0

        for i in range(len(c)):
            sum_c_term += c[i] * (kappa - 2)**(i + 1)

        return np.pi / 12 * (kappa - 2)**4 / (1 + sum_c_term)
    
    elif kappa >= 25:

        sum_d_term = 0
        sum_f_term = 0 

        for i in range(len(d)):
            sum_d_term += d[i] * np.log(kappa)**i 
        
        for i in range(len(f)):
            sum_f_term += f[i] * kappa**-(i + 1)

        return kappa * sum_d_term / (1 - sum_f_term)

# ----------------------------------------------------------------------------------------------------
def phi_02(kappa): # kappa must be an array

    phi = np.zeros_like(kappa)

    mask_1 = kappa < 25
    mask_2 = kappa >= 25

    c1 = 0.8048
    c2 = 0.1459
    c3 = 1.137 * 1.e-3
    c4 = -3.879 * 1.e-6
    c = [c1, c2, c3, c4]

    d0 = -170 + 84*np.log(2) - 16*np.log(2)**2 + np.pi**2/3 * (10 - 4*np.log(2)) + 8*zeta(3)
    d1 = 88 - 40*np.log(2) + 8*np.log(2)**2 - 4/3*np.pi**2
    d2 = -20 + 8*np.log(2)
    d3 = 8/3
    d = [d0, d1, d2, d3]

    f1 = 2.910
    f2 = 78.35
    f3 = 1837
    f = [f1, f2, f3]

    sum_c_term = np.zeros_like(kappa)
    sum_d_term = np.zeros_like(kappa)
    sum_f_term = np.zeros_like(kappa)

    for i in range(len(c)):
        sum_c_term[mask_1] += c[i] * (kappa[mask_1] - 2)**(i + 1)
    
    for i in range(len(d)):
        sum_d_term[mask_2] += d[i] * np.log(kappa[mask_2])**i 
        
    for i in range(len(f)):
        sum_f_term[mask_2] += f[i] * kappa[mask_2]**-(i + 1)

    phi[mask_1] = np.pi / 12 * (kappa[mask_1] - 2)**4 / (1 + sum_c_term[mask_1])
    phi[mask_2] = kappa[mask_2] * sum_d_term[mask_2] / (1 - sum_f_term[mask_2])

    return phi

# ----------------------------------------------------------------------------------------------------
def plot_phi_kappa2():

    kappa = np.logspace(0, 4, num = 100)

    phi_kappa2 = np.zeros_like(kappa)
    for i in range(len(kappa)):
        phi_kappa2[i] = phi_01(kappa[i]) / kappa[i]**2
    
    plt.plot(kappa, phi_kappa2, label = r'LAD ($\kappa$ scalar)') 
    plt.plot(kappa, phi_02(kappa) / kappa**2, label = r'LAD ($\kappa$ array)')

    data = np.loadtxt(f"{REFERENCES_DIR}/Chodorowski1992_Figure2.dat")
    plt.plot(data[:,0], data[:,1], label = 'Chodorowski et al. (1992)')

    plt.xscale('log')
    plt.xlim([1, 1e4])
    plt.ylim([0, 2.5])
    plt.xlabel(r'$\kappa$')
    plt.ylabel(r'$\varphi(\kappa) / \kappa^2$')
    plt.legend()
    plt.savefig(f"{FIGURES_DIR}/Chodorowski1992_xchecks.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/Chodorowski1992_xchecks.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_phi_kappa2()

# ----------------------------------------------------------------------------------------------------