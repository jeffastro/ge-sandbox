import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import cosmology
from scipy import constants
# import warnings; warnings.filterwarnings('ignore')
import sys

def sigma_D(redshift, sigma_z):
    D_H = 3e6/0.7 #kpc
    preamble = (D_H*(redshift)**2)/((1+redshift)**4)
    first = (1+redshift)/np.sqrt(0.3*(1+redshift)**3+0.7)
    second = cosmo.angular_diameter_distance(redshift).value
    return np.sqrt(np.abs(preamble*(first-second)*(sigma_z**2)))
    
    
def phys_radius(redshift, theta, cosmo):
    return cosmo.angular_diameter_distance(redshift).value*theta * 0.0174533 * 1000

def sigma_R(theta, D, sigma_theta, sigma_D):
    return np.sqrt(np.abs((theta**2)*(sigma_D**2)+(D**2)*(sigma_theta**2)+(sigma_D**2)*(sigma_theta**2)))


cosmo = FlatLambdaCDM(H0 =70, Om0 = 0.3)

# ------get data-------------------------------------------------------
thecat = fits.getdata("matched_cs82_redmag_051916.fits")

spher = thecat['SPHEROID_REFF_WORLD'].copy()
redshift = thecat['ZREDMAGIC'].copy()
zband = thecat['ZBAND'].copy()
sigma_zband = thecat['ZBANDERR'].copy()
sigma_ang_r = thecat['SPHEROID_REFFERR_WORLD'].copy()
sigma_z = thecat['ZREDMAGICERR'].copy()


#------compute radii-----------------------------------------------------
radius=[]
'''for gal in xrange(len(spher)):
    grizzly = cosmo.angular_diameter_distance(redshift[gal]).value*spher[gal] * 0.0174533 * 1000
    radius.append(grizzly)'''
for gal in xrange(len(spher)):
    radius.append(phys_radius(redshift[gal], spher[gal], cosmo))
radii = np.array(radius)

print "radii:"
print radii[:10]
print
print
print
print
print
print

#------error propagation-------------------------------------------------
sigmaD = []
'''for gal in xrange(len(spher)):
#    H0=70
    H0 = 2.268308489954634e-18
    D_H = 3e6/0.7 #kpc
    preamble = (D_H*(redshift[gal])**2)/((1+redshift[gal])**4)
    first = (1+redshift[gal])/np.sqrt(0.3*(1+redshift[gal])**3+0.7)
    second = cosmo.angular_diameter_distance(redshift[gal]).value
    sigmaD.append(np.sqrt(np.abs(preamble*(first-second)*(sigma_z[gal]**2))))'''
for gal in xrange(len(spher)):
    sigmaD.append(sigma_D(redshift[gal], sigma_z[gal]))
sigmaD = np.array(sigmaD)

print "sigmaD:"
print sigmaD[:10]
print
print
print
print
print
print

D_A = []
for z in redshift:
    D_A.append(cosmo.angular_diameter_distance(z).value)

D_A = np.array(D_A)

lowz_dist=D_A[(redshift<=(0.26)) & (redshift>=(0.24)) & (radii>1)]                                                                
medz_dist=D_A[(redshift<=(0.51)) & (redshift>=(0.49)) & (radii>1)]  
highz_dist=D_A[(redshift<=(0.67)) & (redshift>=(0.65)) & (radii>1)]  

lowz_err=sigmaD[(redshift<=(0.26)) & (redshift>=(0.24)) & (radii>1)]
medz_err=sigmaD[(redshift<=(0.51)) & (redshift>=(0.49)) & (radii>1)]  
highz_err=sigmaD[(redshift<=(0.67)) & (redshift>=(0.65)) & (radii>1)] 

low_y = lowz_err/lowz_dist
med_y = medz_err/medz_dist
high_y = highz_err/highz_dist
lowrange = xrange(len(low_y))
medrange = np.arange(len(low_y), (len(low_y)+len(med_y)))
highrange = np.arange((len(low_y)+len(med_y)), (len(low_y)+len(med_y)+len(high_y)))

plt.figure()
plt.scatter(lowrange, low_y, color='skyblue', marker=".", label="$z\sim0.25$")
plt.scatter(medrange, med_y, color = 'olive',alpha=0.5, marker=".", label="$z\sim0.5$")
plt.scatter(highrange, high_y, color = 'salmon',alpha=0.5, marker=".", label="$z\sim0.66$")
plt.ylabel("$\sigma_{D_A}/D_A$")
plt.yscale('log')
plt.legend(loc="upper left")
plt.title("Fractional error in angular diameter distance")
# plt.axis([0, len(low_y)+len(med_y),0,.25])
plt.show()


sigmaR = []
for gal in xrange(len(spher)):
    sigmaR.append(sigma_R(spher[gal], cosmo.angular_diameter_distance(redshift[gal]).value, sigma_ang_r[gal], sigmaD[gal]))
sigmaR=np.array(sigmaR)
print "sigmaR:"
print sigmaR[:10]
print
print
print
print
print
print
    
#-----05.19.16 stoppped here-----


#ADD HORIZONTAL ERROR BARS

lowz_radii=radii[(redshift<=(0.26)) & (redshift>=(0.24)) & (radii>1)]                                                                
medz_radii=radii[(redshift<=(0.51)) & (redshift>=(0.49)) & (radii>1)]  
highz_radii=radii[(redshift<=(0.67)) & (redshift>=(0.65)) & (radii>1)]  

lowz_err=sigmaR[(redshift<=(0.26)) & (redshift>=(0.24)) & (radii>1)]
medz_err=sigmaR[(redshift<=(0.51)) & (redshift>=(0.49)) & (radii>1)]  
highz_err=sigmaR[(redshift<=(0.67)) & (redshift>=(0.65)) & (radii>1)] 

low_y = lowz_err/lowz_radii
med_y = medz_err/medz_radii
high_y = highz_err/highz_radii
lowrange = xrange(len(low_y))
medrange = np.arange(len(low_y), (len(low_y)+len(med_y)))
highrange = np.arange((len(low_y)+len(med_y)), (len(low_y)+len(med_y)+len(high_y)))


plt.figure()
plt.scatter(lowrange, low_y, color='b', marker=".", label="$z\sim0.25$")
plt.scatter(medrange, med_y, color = 'g',alpha=0.5, marker=".", label="$z\sim0.5$")
plt.scatter(highrange, high_y, color = 'magenta',alpha=0.5, marker=".", label="$z\sim0.66$")
plt.ylabel("$\sigma_R/radius$")
plt.yscale('log')
plt.legend(loc="upper left")
plt.title("Fractional error in physical radius")
plt.savefig("sigmaR_over_R.png")
# plt.axis([0, len(low_y)+len(med_y),0,.25])
plt.show()
