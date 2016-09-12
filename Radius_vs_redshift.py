import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM


cosmo = FlatLambdaCDM(H0 =70, Om0 = 0.3)


thecat = fits.getdata('matched_cs82_redmag_040316.fits')
emceedata = fits.getdata("emceeRun031616_quad.fits")

intercepts = emceedata['a'].copy()
lins = emceedata['b'].copy()
quads = emceedata['c'].copy()

sigmas = emceedata['sigma'].copy()
factors = emceedata['factor'].copy()
emceeredshifts = emceedata['redshift'].copy()
galnumbers = emceedata['number_of_galaxies'].copy()

spher = thecat['SPHEROID_REFF_WORLD'].copy()
redshift = thecat['ZREDMAGIC'].copy()
zband = thecat['ZBAND'].copy()

#get pivot point
kavli_fits = fits.getdata("matchedkavli112315.fits")
zbandforpivot = kavli_fits['ZBAND'].copy()
pivot = np.median(zbandforpivot)
print pivot
# np.savetxt("zband_pivot_matchedkavli112315.txt", np.array(pivot))

radius = []
mass = np.zeros(len(zband)).tolist()
# print mass

width = 0.02
halfwidth = width / 2

#convert angular size to physical size
for gal in xrange(len(spher)):
    grizzly = cosmo.angular_diameter_distance(redshift[gal]).value*spher[gal] * 0.0174533 * 1000
    radius.append(grizzly)

for gal in xrange(len(zband)):
	currentz = redshift[gal]
	mcz = min(emceeredshifts, key=lambda x:abs(x-currentz))
	mass[gal]=intercepts[emceeredshifts==mcz] + lins[emceeredshifts==mcz]*(zband[gal]-pivot) + quads[emceeredshifts==mcz]*(zband[gal]-pivot)**2

mass=np.array(mass).flatten()
print mass
print mass.shape

radius=np.array(radius)


radii_10_8=[]
radii_11_1=[]
radii_11_4=[]

currentmstar = 10.8
for zbin in emceeredshifts:
	currentradius = radius[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20) ].flatten()
	massofinterest = mass[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20)].flatten()
	currentfit = np.polyfit(massofinterest, currentradius, 2)
	currentquad,currentlin,currentinter=np.poly1d(currentfit)
	bestfitradius = currentinter + currentlin * (currentmstar ) + currentquad * (currentmstar) ** 2
	radii_10_8.append(bestfitradius)

currentmstar = 11.1
for zbin in emceeredshifts:
	currentradius = radius[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20)].flatten()
	massofinterest = mass[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20)].flatten()
	currentfit = np.polyfit(massofinterest, currentradius, 2)
	currentquad,currentlin,currentinter=np.poly1d(currentfit)
	bestfitradius = currentinter + currentlin * (currentmstar ) + currentquad * (currentmstar) ** 2
	radii_11_1.append(bestfitradius)
	
currentmstar = 11.4
for zbin in emceeredshifts:
	currentradius = radius[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20) ].flatten()
	massofinterest = mass[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20)].flatten()
	currentfit = np.polyfit(massofinterest, currentradius, 2)
	currentquad,currentlin,currentinter=np.poly1d(currentfit)
	bestfitradius = currentinter + currentlin * (currentmstar ) + currentquad * (currentmstar) ** 2
	radii_11_4.append(bestfitradius)

lowfit = np.polyfit(emceeredshifts, radii_10_8, 2)
lowquad , lowlin, lowinter = (lowfit)
bestlowfit = lowinter + lowlin * (emceeredshifts ) + lowquad * (emceeredshifts) ** 2

midfit = np.polyfit(emceeredshifts, radii_11_1, 2)
midquad , midlin, midinter = (midfit)
bestmidfit = midinter + midlin * (emceeredshifts ) + midquad * (emceeredshifts) ** 2

highfit = np.polyfit(emceeredshifts, radii_11_4, 2)
highquad , highlin, highinter = (highfit)
besthighfit = highinter + highlin * (emceeredshifts ) + highquad * (emceeredshifts) ** 2

	
plt.figure(facecolor='white')
plt.scatter(emceeredshifts, radii_10_8, marker='o', edgecolor='none',color='darkseagreen', label="$10.8$ $Log$ $M_\odot$")
plt.scatter(emceeredshifts, radii_11_1, marker='o', edgecolor='none',color='dodgerblue', label="$11.1$ $Log$ $M_\odot$")
plt.scatter(emceeredshifts, radii_11_4, marker='o', edgecolor='none',color='darkorchid', label="$11.4$ $Log$ $M_\odot$")
plt.plot(emceeredshifts, bestlowfit,color='darkseagreen')
plt.plot(emceeredshifts, bestmidfit,color='dodgerblue')
plt.plot(emceeredshifts, besthighfit,color='darkorchid')
plt.legend(loc="lower left")
plt.ylabel("$Physical$ $radius$ $(kpc)$")
plt.xlabel("$Redshift$")
plt.title("$Radius$ $vs$ $redshift$ $in$ $stellar$ $mass$ $bins$")
plt.xlim([0.2,0.7])
plt.show()