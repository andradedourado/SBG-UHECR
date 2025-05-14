from astropy.cosmology import Planck18 as cosmo
from astropy.cosmology import z_at_value
from astropy import units as u
import numpy as np

RESULTS_DIR = "../results"

data = np.genfromtxt(f'starburst_galaxies.dat', dtype = None, encoding = None)

# ----------------------------------------------------------------------------------------------------
def lum_dist_to_z():

    d_lum, z = np.zeros(len(data)), np.zeros(len(data))
    
    for idata in range(len(data)):
        d_lum[idata] = data[idata][8]
        z[idata] = z_at_value(cosmo.luminosity_distance, d_lum[idata] * u.Mpc) 

    return d_lum, z

# ----------------------------------------------------------------------------------------------------
def save_lum_dist_to_z():

    d_lum, z = lum_dist_to_z()
    np.savetxt(f"{RESULTS_DIR}/lum_dist_to_z.dat", np.column_stack((d_lum, z)), fmt = "%.15e")
 
# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    save_lum_dist_to_z()

# ----------------------------------------------------------------------------------------------------