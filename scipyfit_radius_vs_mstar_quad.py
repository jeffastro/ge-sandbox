# Last modified 03.27.16 at 19:35
# loop for creating many quadratic fits across z bins with a single pivot for all.
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import corner
import emcee
from astropy.cosmology import FlatLambdaCDM
import warnings; warnings.filterwarnings('ignore')
# prettyplotlib imports 
import prettyplotlib as ppl
import matplotlib as mpl
from prettyplotlib import brewer2mpl


cosmo = FlatLambdaCDM(H0 =70, Om0 = 0.3)

# ------get data-------------------------------------------------------
thecat = fits.getdata("matched_cs82_redmag_mgc_031816.fits")
mstar = thecat['MASS_IR_BEST'].copy()
spher = thecat['SPHEROID_REFF_WORLD'].copy()
redshift = thecat['ZREDMAGIC'].copy()
radius = []

#convert angular size to physical size
for gal in xrange(len(spher)):
    grizzly = cosmo.angular_diameter_distance(redshift[gal]).value*spher[gal] * 0.0174533 * 1000
    radius.append(grizzly)

radius = np.array(radius)
zbin = 0.25
maxredshift=0.65
width = 0.1
halfwidth = width / 2
pivot = np.median(mstar)

n = 0

a = []
alow = []
ahigh = []
b = []
blow = []
bhigh = []
c = []
clow = []
chigh = []
sigmas = []
siglow = []
sighigh = []
factors = []
factorlow = []
factorhigh = []
redshifts = []
galnumbers = []

while (zbin < maxredshift):
    medz_mstar = mstar[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (mstar > 0)  & (radius>0)]
    medz_radius = radius[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (mstar > 0) & (radius>0)]
    fit = np.polyfit(medz_mstar, medz_radius, 2)
    quad,lin,inter=np.poly1d(fit)
    a.append(inter)
    b.append(lin)
    c.append(quad)
    redshifts.append(zbin)
    n+=1
    zbin+=width



lines = zip(a,b,c, redshifts)

xranges = np.arange(10.6, 11.4, 0.1)
fig, ax = plt.subplots(1)
for fit in lines:
	print fit
	a, b, c, z = fit
	bestfit = a + b * (xranges ) + c * (xranges) ** 2
	ppl.plot(ax, xranges, bestfit, label=str(z), linewidth=0.75)

ppl.legend(ax, loc='upper left', ncol=3)

plt.title("Galaxy radius vs mass")
plt.xlabel("Log solar mass")
plt.ylabel("Physical radius (kpc)")
ax.set_yscale('log')
plt.show()