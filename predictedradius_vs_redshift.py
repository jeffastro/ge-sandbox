#Last edited 03.16.16 at 00:45
#plots predicted Radius(at ZBAND=19.5) vs redshift                                                               

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import pandas as pd
import sys
from matplotlib import cm

catalog = fits.getdata("matched_cs82_redmag_mgc_031816.fits")
emcee = fits.getdata("emceeRun022816_linear.fits")
cosmo = FlatLambdaCDM(H0 =70, Om0 = 0.3)

mstar = catalog['MASS_IR_BEST'].copy()
spheroid = catalog['SPHEROID_REFF_WORLD'].copy()
zband = catalog['ZBAND'].copy()
redshift=catalog['ZREDMAGIC'].copy()
mcredshifts=emcee['redshift'].copy()
slopes = emcee['a'].copy()
intercepts=emcee['b'].copy()
radii = []

for gal in xrange(len(spheroid)):
    grizzly = cosmo.angular_diameter_distance(redshift[gal]).value*spheroid[gal] * 0.0174533 * 1000
    radii.append(grizzly)

averageredshifts = []
predictedradii = []
predictedmass = []
pivot = np.median(zband)
radii=np.array(radii)
n=7
print "emcee redshifts:"
print mcredshifts
mcredshifts = mcredshifts[n:]
centralzband = 19

for zused in mcredshifts:
	zbandofinterest=zband[(redshift<=(zused+0.02)) & (redshift>=(zused)) ]

	radiiofinterest=radii[(redshift<=(zused+0.02)) & (redshift>=(zused))]
	slope = slopes[n]
	intercept = intercepts[n]
	if (len(zbandofinterest)>0):
		zbandofinterest=np.array(zbandofinterest)
		radiiofinterest=np.array(radiiofinterest)
		fit = np.polyfit(zbandofinterest, radiiofinterest, 1)
		coucou=np.poly1d(fit)
		pradius=coucou(centralzband) 
		pmass=slope * (centralzband-pivot) + intercept
		predictedradii.append(pradius)
		predictedmass.append(pmass)
	else:
		print "zband has nothing"
	n+=1

print predictedmass
print predictedradii
plt.figure()
plt.scatter(predictedmass,predictedradii,c=mcredshifts, cmap=cm.get_cmap('rainbow'))
plt.title ("Radius (cs82) at ZBAND=%d (redMaGiC) vs Stellar Mass (MGC)\nPredicted using linear emcee fit" % centralzband)
plt.ylabel("Predicted Radius (kpc) at ZBAND=%d" % centralzband)
plt.xlabel("Predicted M* (MGC)")
plt.colorbar()
plt.yscale('log', basey=10)
plt.show()