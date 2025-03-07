import matplotlib.pyplot as plt
import numpy as np 

# ----------------------------------------------------------------------------------------------------
def inelasticity(eps_prime):

    Y_inf = 0.47
    xb = 6e9
    dlt = 0.33
    s = 0.15

    return 

# ----------------------------------------------------------------------------------------------------
def plot_photopion_inelasticity():

    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ is '__main__':

    plot_photopion_inelasticity()

# ----------------------------------------------------------------------------------------------------




# double inelasticity(double epsPrime) {
#   static const double Y0 = 0.47;
#   static const double b = 6e9 * SI::eV;
#   static const double d = 0.33;
#   static const double s = 0.15;
#   const auto x = epsPrime / b;
#   return Y0 * std::pow(x, d) / std::pow(1. + std::pow(x, d / s), s);
# }