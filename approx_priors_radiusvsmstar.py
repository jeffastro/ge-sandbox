#Last updated 3.27.15 at 15:44

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import warnings; warnings.filterwarnings('ignore')

catalog = fits.getdata("kavli_cs82_113015.fits")
cosmo = FlatLambdaCDM(H0 =70, Om0 = 0.3)

mstar = catalog['MASS_IR_BEST'].copy()
zbest = catalog['ZBEST'].copy()
spheroid = catalog['SPHEROID_REFF_WORLD'].copy()
radius = []
pivot = np.median(mstar)

for gal in xrange(len(spheroid)):
    grizzly = cosmo.angular_diameter_distance(zbest[gal]).value*spheroid[gal] * 0.0174533 * 1000
    radius.append(grizzly)

radius = np.array(radius)
# print radius
# lowzrad = radius[(zbest<0.251) & (zbest>0.249)]
medzrad = radius[(zbest<0.501) & (zbest>0.499) & (radius>0)]-pivot
# highzrad = radius[(zbest<0.701) & (zbest>0.699)]
# lowmstar = mstar[(zbest<0.251) & (zbest>0.249)]
medmstar = mstar[(zbest<0.501) & (zbest>0.499) & (radius>0)]
# highmstar = mstar[(zbest<0.701) & (zbest>0.699)]

fit = np.polyfit(medmstar, medzrad, 2)
coucou=np.poly1d(fit)
print coucou

medfit=coucou(medmstar)
# print medfit

plt.figure()
plt.plot(medmstar, medfit)
# plt.scatter(lowmstar, lowzrad, label="z~0.25", color = 'm', alpha = 0.2)
# plt.scatter(medmstar, medzrad, label="z~0.50", color = 'g', alpha = 0.1)
# plt.scatter(highmstar, highzrad, label="z~0.70", color = 'b', alpha = 0.05)
# plt.axis([9,12.5,1,300])
# plt.title(# "Galaxy radius (cs82) vs M* (Massive Galaxy Catalog)")
# plt.xlabel("MASS_IR_BEST (Log of a thousand solar masses)")
# plt.ylabel("Radius (kpc)")
# # plt.yscale('log')
# plt.legend(loc='best')
plt.show()
# plt.savefig('radius_vs_mstar_redshift_bins.pdf')