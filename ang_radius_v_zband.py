import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
# import warnings; warnings.filterwarnings('ignore')

cosmo = FlatLambdaCDM(H0 =70, Om0 = 0.3)

# ------get data-------------------------------------------------------
thecat = fits.getdata("matched_cs82_redmag_040316.fits")

spher = thecat['SPHEROID_REFF_WORLD'].copy()
redshift = thecat['ZREDMAGIC'].copy()
zband = thecat['ZBAND'].copy()

radius=[]

'''for gal in xrange(len(spher)):
    grizzly = cosmo.angular_diameter_distance(redshift[gal]).value*spher[gal] * 0.0174533 * 1000
    radius.append(grizzly)

radii = np.array(radius)'''

lowz_radii=spher[(redshift<=(0.26)) & (redshift>=(0.24)) & (spher>1e-6)]                                                                
medz_radii=spher[(redshift<=(0.51)) & (redshift>=(0.49))& (spher>1e-6)]  
highz_radii=spher[(redshift<=(0.67)) & (redshift>=(0.65))& (spher>1e-6)]  

lowz_zband=zband[(redshift<=(0.26)) & (redshift>=(0.24))& (spher>1e-6)]                                                                
medz_zband=zband[(redshift<=(0.51)) & (redshift>=(0.49))& (spher>1e-6)]  
highz_zband=zband[(redshift<=(0.67)) & (redshift>=(0.65))& (spher>1e-6)]  
                                                                                      
pictures=plt.figure(facecolor='white')
pic=plt.plot(markersize=3)
#lowzfunc = stats.linregress(lowz[:,2], lowz[:,3])

#these are the prepackaged fit functions
lowzfunc = np.polyfit(lowz_zband, lowz_radii, 1)
medzfunc = np.polyfit(medz_zband, medz_radii, 1)
highzfunc = np.polyfit(highz_zband, highz_radii, 1)
x=np.linspace (zband.min(), zband.max(), 100)


p1 = np.poly1d(lowzfunc)
p2 = np.poly1d(medzfunc)
p3 = np.poly1d(highzfunc)
plt.plot(x, p1(x), color='b')
plt.plot(x, p2(x), color='g')
plt.plot(x, p3(x), color='m')


plt.scatter(lowz_zband, lowz_radii, color='b', marker="*", alpha=0.15, label = '$z\ 0.25$')
plt.scatter(medz_zband, medz_radii, color='g', marker="*", alpha=0.13, label = '$z\ 0.5$')
plt.scatter(highz_zband, highz_radii, color='m', marker="*", alpha=0.1,label = "$z\ 0.66$")


plt.title ("$Angular$ $radius$ $(cs82)$ $vs$ $apparent$ $magnitude$ $in$ $z$ $band$ $(stripe82)$")
plt.xlabel("$Apparent$ $magnitude$ $in$ $z$ $band$")
plt.ylabel("$Angular$ $Radius$")
#plt.xscale('log', basex=10)
plt.yscale('log', basey=10)
plt.legend(loc="lower left")

plt.axis([15.8, 22.5,5e-5,1e-2])
plt.show()