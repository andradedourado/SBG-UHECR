from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'legend.fontsize': 'medium',
'legend.title_fontsize': 'medium',
'axes.labelsize': 'large',
'axes.titlesize': 'x-large',
'xtick.labelsize': 'large',
'ytick.labelsize': 'large'})

DATA_DIR = "../data"
FIGURES_DIR = "../figures"

data = np.genfromtxt(f'starburst_galaxies.dat', dtype = None, encoding = None)
event_data = np.loadtxt(f'{DATA_DIR}/AugerApJS2022_Yr_JD_UTC_Th_Ph_RA_Dec_E_Expo.dat')

degree = np.pi / 180.

# ----------------------------------------------------------------------------------------------------
def countours_auger_sky(nbins = int(1.e2)):
	
    dec = np.full(nbins, 45) 
    ra = np.linspace(0, 360, num = nbins) 

    countours = np.zeros((nbins, 2))
    
    for icoord in range(nbins):

        c_equatorial = SkyCoord(ra = ra[icoord] * u.degree, dec = dec[icoord] * u.degree, frame = 'icrs')
        c_galactic = c_equatorial.galactic
		
        l = c_galactic.l.radian
        b = c_galactic.b.radian
		
        # Python coordinates
        if 2*np.pi >= l >= np.pi:
            l = 2*np.pi - l
        elif 0 <= l <= np.pi:
            l = -l
        
        countours[icoord, 0] = b
        countours[icoord, 1] = l

    return countours

# ----------------------------------------------------------------------------------------------------
def ra_dec_to_gal_lat(ra, dec):
     
    c_equatorial = SkyCoord(ra = ra * u.degree, dec = dec * u.degree, frame = 'icrs', unit = 'deg')
    c_galactic = c_equatorial.galactic
			
    return c_galactic.b.radian

# ----------------------------------------------------------------------------------------------------
def ra_dec_to_gal_lon(ra, dec):

    c_equatorial = SkyCoord(ra = ra * u.degree, dec = dec * u.degree, frame = 'icrs', unit = 'deg')
    c_galactic = c_equatorial.galactic
			
    l = c_galactic.l.radian	

    l_transformed = np.empty_like(l)

    # Python coordinates
    mask1 = (l >= np.pi) & (l <= 2 * np.pi)
    mask2 = (l >= 0) & (l <= np.pi)

    l_transformed[mask1] = 2 * np.pi - l[mask1]
    l_transformed[mask2] = -l[mask2]

    return l_transformed		

# ----------------------------------------------------------------------------------------------------
def get_directions_in_galactic_coords(): # Galaxies

    gal_lat = np.zeros(len(data))
    gal_lon = np.zeros(len(data))
    dist_mpc = np.zeros(len(data))

    for idata in range(len(data)):
        gal_lat[idata] = ra_dec_to_gal_lat(data[idata][4], data[idata][5])
        gal_lon[idata] = ra_dec_to_gal_lon(data[idata][4], data[idata][5])
        dist_mpc[idata] = data[idata][8]

    return gal_lat, gal_lon, dist_mpc

# ----------------------------------------------------------------------------------------------------
def plot_galaxy_directions():

    plt.figure()
    axs = plt.subplot(111, projection = 'mollweide')
    axs.set_longitude_grid_ends(90)

    plt.fill(countours_auger_sky()[:,1], countours_auger_sky()[:,0], facecolor = 'lightgray', edgecolor = None)

    gal_lat, gal_lon, dist_mpc = get_directions_in_galactic_coords()
    s_factor = 1000 
    marker_sizes = s_factor / (dist_mpc ** 2)
    plt.scatter(gal_lon, gal_lat, s = marker_sizes, edgecolors = 'k', linewidths = 0.675, c = 'orangered')
    
    plt.xlabel(r'Galactic longitude, $l \: {\rm [deg]}$', labelpad = 10)
    plt.ylabel(r'Galactic latitude, $b \: {\rm [deg]}$')
    plt.xticks(ticks = [-120 * degree, -60 * degree, 0 * degree, 60 * degree, 120 * degree],
    labels = [r'$120\degree$', r'$60\degree$', r'$0\degree$', r'$300\degree$', r'$240\degree$'])
    plt.yticks(ticks = [-60 * degree, -30 * degree, 0 * degree, 30 * degree, 60 * degree], 
	labels = [r'$-60\degree$', r'$-30\degree$', r'$0\degree$', r'$30\degree$', r'$60\degree$'])
    plt.grid(linestyle = 'dotted', color = 'black', linewidth = 0.5, zorder = -1.0)
    plt.savefig(f'{FIGURES_DIR}/starburst_galaxies.pdf', bbox_inches = 'tight')
    plt.savefig(f'{FIGURES_DIR}/starburst_galaxies.png', bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
def plot_galaxy_directions_wEvents():

    plt.figure()
    axs = plt.subplot(111, projection = 'mollweide')
    axs.set_longitude_grid_ends(90)

    plt.fill(countours_auger_sky()[:,1], countours_auger_sky()[:,0], facecolor = 'lightgray', edgecolor = None)

    plt.scatter(ra_dec_to_gal_lon(event_data[:,5], event_data[:,6]), ra_dec_to_gal_lat(event_data[:,5], event_data[:,6]), alpha = 0.05, c = 'gray', marker = '.') # edgecolors = 'none'

    gal_lat, gal_lon, dist_mpc = get_directions_in_galactic_coords()
    s_factor = 1000 
    marker_sizes = s_factor / (dist_mpc ** 2)
    plt.scatter(gal_lon, gal_lat, s = marker_sizes, edgecolors = 'k', linewidths = 0.675, c = 'orangered')

    plt.xlabel(r'Galactic longitude, $l \: {\rm [deg]}$', labelpad = 10)
    plt.ylabel(r'Galactic latitude, $b \: {\rm [deg]}$')
    plt.xticks(ticks = [-120 * degree, -60 * degree, 0 * degree, 60 * degree, 120 * degree],
    labels = [r'$120\degree$', r'$60\degree$', r'$0\degree$', r'$300\degree$', r'$240\degree$'])
    plt.yticks(ticks = [-60 * degree, -30 * degree, 0 * degree, 30 * degree, 60 * degree], 
	labels = [r'$-60\degree$', r'$-30\degree$', r'$0\degree$', r'$30\degree$', r'$60\degree$'])
    plt.grid(linestyle = 'dotted', color = 'black', linewidth = 0.5, zorder = -1.0)
    plt.savefig(f'{FIGURES_DIR}/starburst_galaxies_wEvents.pdf', bbox_inches = 'tight')
    plt.savefig(f'{FIGURES_DIR}/starburst_galaxies_wEvents.png', bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    # plot_galaxy_directions()
    plot_galaxy_directions_wEvents()

# ----------------------------------------------------------------------------------------------------