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

data = np.genfromtxt('starburst_galaxies.dat', dtype = None, encoding = None)
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
def ra_dec_to_gal_lat(idata):
     
    c_equatorial = SkyCoord(ra = data[idata][4] * u.degree, dec = data[idata][5] * u.degree, frame = 'icrs', unit = 'deg')
    c_galactic = c_equatorial.galactic
			
    return c_galactic.b.radian

# ----------------------------------------------------------------------------------------------------
def ra_dec_to_gal_lon(idata):

    c_equatorial = SkyCoord(ra = data[idata][4] * u.degree, dec = data[idata][5] * u.degree, frame = 'icrs', unit = 'deg')
    c_galactic = c_equatorial.galactic
			
    l = c_galactic.l.radian			
			
    if 2*np.pi >= l >= np.pi:			
        return 2*np.pi - l
    elif 0 <= l <= np.pi:
	    return -l

# ----------------------------------------------------------------------------------------------------
def plot_galaxy_directions():

    plt.figure()
    axs = plt.subplot(111, projection = 'mollweide')
    axs.set_longitude_grid_ends(90)

    gal_lat = np.zeros(len(data))
    gal_lon = np.zeros(len(data))
    dist_mpc = np.zeros(len(data))

    for idata in range(len(data)):
        gal_lat[idata] = ra_dec_to_gal_lat(idata)
        gal_lon[idata] = ra_dec_to_gal_lon(idata)
        dist_mpc[idata] = data[idata][8]

    s_factor = 1000 
    marker_sizes = s_factor / (dist_mpc ** 2)
	
    plt.fill(countours_auger_sky()[:,1], countours_auger_sky()[:,0], facecolor = 'lightgray', edgecolor = None)
    plt.scatter(gal_lon, gal_lat, s = marker_sizes, edgecolors = 'k', linewidths = 0.675, c = 'orangered')
    
    plt.xlabel(r'Galactic longitude, $l \: {\rm [deg]}$', labelpad = 10)
    plt.ylabel(r'Galactic latitude, $b \: {\rm [deg]}$')
    plt.xticks(ticks = [-120 * degree, -60 * degree, 0 * degree, 60 * degree, 120 * degree],
    labels = [r'$120\degree$', r'$60\degree$', r'$0\degree$', r'$300\degree$', r'$240\degree$'])
    plt.yticks(ticks = [-60 * degree, -30 * degree, 0 * degree, 30 * degree, 60 * degree], 
	labels = [r'$-60\degree$', r'$-30\degree$', r'$0\degree$', r'$30\degree$', r'$60\degree$'])
    plt.grid(linestyle = 'dotted', color = 'black', linewidth = 0.5, zorder = -1.0)
    plt.savefig('starburst_galaxies.pdf', format = 'pdf', bbox_inches = 'tight')
    plt.savefig('starburst_galaxies.png', format = 'png', bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_galaxy_directions()

# ----------------------------------------------------------------------------------------------------

