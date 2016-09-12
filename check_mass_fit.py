import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import warnings; warnings.filterwarnings('ignore')


# ------get data-------------------------------------------------------
thecat = fits.open("matched_cs82_redmag_mgc_031716.fits")
# emceedata = fits.getdata("emceeRun031616_quad.fits")
emceedata = fits.open("emceeRun022816_linear.fits")

intercepts = emceedata['b'].copy()
lins = emceedata['a'].copy()
# intercepts = emceedata['a'].copy()
# lins = emceedata['b'].copy()
# quads = emceedata['c'].copy()
emceeredshifts = emceedata['redshift'].copy()

mstar_kavli = thecat['MASS_IR_BEST'].copy()
redshift = thecat['ZREDMAGIC'].copy()
zband = thecat['ZBAND'].copy()
#print max(redshift)	
thecat.close()
emceedata.close()

#get pivot point
kavli_fits = fits.open("matchedkavli112315.fits")
zbandforpivot = kavli_fits['ZBAND'].copy()
pivot = np.median(zbandforpivot)
# print pivot
kavli_fits.close()

mass = np.empty(len(zband)).tolist()

width = 0.02
halfwidth = width / 2

#---test unit----
# gal=10
# currentz = redshift[gal]
# print "currentz: " + str(currentz)
# mcz = min(emceeredshifts, key=lambda x:abs(x-currentz))
# print "mcz: " + str(mcz)
# print "intercept: " + str(intercepts[emceeredshifts==mcz])
# print "intercept by hand: "  + str(intercepts[16])
# print "lin: " + str(lins[emceeredshifts==mcz])
# print "lin by hand: " + str(lins[16])
# print "automatic lin mult: " +str(lins[emceeredshifts==mcz]*(zband[gal]-pivot))
# print "lin mult by hand: " + str(lins[16]*(zband[gal]-pivot))
# mass[gal]=intercepts[emceeredshifts==mcz] + lins[emceeredshifts==mcz]*(zband[gal]-pivot) + quads[emceeredshifts==mcz]*(zband[gal]-pivot)**2
# print mass[gal]

for gal in xrange(10):#(len(zband)):
	currentz = redshift[gal]
	mcz = min(emceeredshifts, key=lambda x:abs(x-currentz))
# 	print intercepts[emceeredshifts==mcz]
	mass[gal]=intercepts[emceeredshifts==mcz] + lins[emceeredshifts==mcz]*(zband[gal]-pivot) #+ quads[emceeredshifts==mcz]*(zband[gal]-pivot)**2

mass=np.array(mass)
newmass = mass[(mstar_kavli>0) & (mass>0)]
newmstar_kavli = mstar_kavli[(mstar_kavli>0) & (mass>0)]
plt.figure()
plt.scatter(newmstar_kavli[:10], newmass[:10])

plt.show()