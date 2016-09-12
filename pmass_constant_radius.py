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
predictedzband = []
predictedmass = []
pivot = np.median(zband)
radii=np.array(radii)
n=7
print "emcee redshifts:"
print mcredshifts
mcredshifts = mcredshifts[n:]
#centralzband = 19

for zused in mcredshifts:
	zbandofinterest=zband[(redshift<=(zused+0.02)) & (redshift>=(zused)) ]

	radiiofinterest=radii[(redshift<=(zused+0.02)) & (redshift>=(zused))]
# 	centralzband = np.median(zbandofinterest)
	slope = slopes[n]
	intercept = intercepts[n]
	if (len(zbandofinterest)>0):
		zbandofinterest=np.array(zbandofinterest)
		radiiofinterest=np.array(radiiofinterest)
		fit = np.polyfit(radiiofinterest,zbandofinterest, 1)
		coucou=np.poly1d(fit)
		pzband=coucou(5) 
		pmass=slope * (pzband-pivot) + intercept
		predictedzband.append(pzband)
		predictedmass.append(pmass)
	else:
		print "zband has nothing"
	n+=1

print predictedmass
print predictedzband
plt.figure()
plt.scatter(predictedzband,predictedmass,c=mcredshifts, cmap=cm.get_cmap('rainbow'))
plt.title ("Stellar Mass (MGC) vs ZBAND (redMaGiC) predicted at 5kpc radius (cs82)\nlinear emcee fit" )
plt.ylabel("Predicted M* (from MGC)")
plt.xlabel("Predicted ZBAND")
plt.colorbar(label="redshift")
# plt.yscale('log', basey=10)
plt.show()