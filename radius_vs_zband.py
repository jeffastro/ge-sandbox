import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import cosmology
from scipy import constants
# import warnings; warnings.filterwarnings('ignore')

def sigma_D(redshift, sigma_z):
#     H0=70
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

print "H_0 at 0.2: " + cosmology.H(0.2)
print "H_0 value at 0.0: " + cosmology.H(0).value
print "H_0 units at 0.0: " + cosmology.H0(0).unit

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
    H0=70
    D_H = 3e6/H0 #Mpc
    preamble = (D_H*(redshift[gal])**2)/((1+redshift[gal])**4)
    first = (1+redshift[gal])/np.sqrt(0.3*(1+redshift[gal])**3+0.7)
    second = cosmo.angular_diameter_distance(redshift[gal]).value
    sigmaD.append(np.sqrt(np.abs(preamble*(first-second)*(sigma_z[gal]**2))))
'''
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


lowz_zband=zband[(redshift<=(0.26)) & (redshift>=(0.24)) & (radii>1)]                                                                
medz_zband=zband[(redshift<=(0.51)) & (redshift>=(0.49)) & (radii>1)]  
highz_zband=zband[(redshift<=(0.67)) & (redshift>=(0.65)) & (radii>1)]  

lowz_err=sigmaR[(redshift<=(0.26)) & (redshift>=(0.24)) & (radii>1)]
medz_err=sigmaR[(redshift<=(0.51)) & (redshift>=(0.49)) & (radii>1)]  
highz_err=sigmaR[(redshift<=(0.67)) & (redshift>=(0.65)) & (radii>1)] 

lowz_zband_err=sigma_zband[(redshift<=(0.26)) & (redshift>=(0.24)) & (radii>1)]
print lowz_zband_err
medz_zband_err=sigma_zband[(redshift<=(0.51)) & (redshift>=(0.49)) & (radii>1)]  
highz_zband_err=sigma_zband[(redshift<=(0.67)) & (redshift>=(0.65)) & (radii>1)] 
                                                                                      
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


# plt.scatter(lowz_zband, lowz_radii, color='b', marker="*", alpha=0.15, label = '$z\ 0.25$')
# plt.errorbar(lowz_zband, lowz_ radii, yerr = lowz_err,   color='b', capsize=0, alpha=0.8,linestyle="None", label = '$z~0.25$')
# plt.errorbar(medz_zband, medz_radii,yerr = medz_err,  color='g',capsize=0, alpha=0.2,linestyle="None", label = '$z~0.5$')
# plt.errorbar(highz_zband, highz_radii,yerr = highz_err, color='magenta', capsize=0, alpha=0.2,linestyle="None",label = "$z~0.66$")



plt.errorbar(lowz_zband, lowz_radii, yerr = lowz_err, xerr=lowz_zband_err, color='b',capsize=0, alpha=1,linestyle="None", label = '$z\sim0.25$')
plt.errorbar(medz_zband, medz_radii,yerr = medz_err, xerr=medz_zband_err, color='g',capsize=0, alpha=0.5,linestyle="None", label = '$z\sim0.5$')
plt.errorbar(highz_zband, highz_radii,yerr = highz_err, xerr=highz_zband_err, color='magenta',capsize=0, alpha=0.5,linestyle="None",label = "$z\sim0.66$")


plt.title ("$Galaxy$ $radius$ $(cs82)$ $vs$ $apparent$ $magnitude$ $in$ $z$ $band$ $(stripe82):$\n$De$ $Vaucouleurs$ $profile$")
plt.xlabel("$Apparent$ $magnitude$ $in$ $z$ $band$")
plt.ylabel("$Radius$ $(kpc)$")
#plt.xscale('log', basex=10)
plt.yscale('log', basey=10)
plt.legend(loc="lower left")

plt.axis([15.8, 22.5,1,80])
plt.show()
