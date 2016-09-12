#Modified 4.3.16 at 16:17

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import warnings; warnings.filterwarnings('ignore')

cosmo = FlatLambdaCDM(H0 =70, Om0 = 0.3)

# ------get data-------------------------------------------------------
thecat = fits.getdata("matched_cs82_redmag_040316.fits")
# emceedata = fits.getdata("emceeRun022816_linear.fits")
emceedata = fits.getdata("emceeRun031616_quad.fits")

intercepts = emceedata['a'].copy()
lins = emceedata['b'].copy()
quads = emceedata['c'].copy()
# intercepts = emceedata['b'].copy()
# lins = emceedata['a'].copy()

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
# print mass
# print mass.shape

radius=np.array(radius)
# print radius
# print radius.shape
# print mass

# oldhalfwidth=halfwidth
# halfwith=0.001

oldgrayscalecolors=['#D9D9D9',
'#D6D6D6',
'#D4D4D4',
'#CCCCCC',
'#C9C9C9',
'#C7C7C7',
'#BFBFBF',
'#BDBDBD',
'#B5B5B5',
'#B3B3B3',
'#ABABAB',
'#A8A8A8',
'#A6A6A6',
'#A3A3A3',
'#9C9C9C',
'#999999',
'#969696',
'#949494',
'#919191',
'#8F8F8F',
'#8C8C8C',
'#8A8A8A',
'#878787',
'#858585',
'#828282',
'#7F7F7F',
'#7D7D7D',
'#7A7A7A',
'#787878',
'#757575',
'#737373',
'#707070',
'#6E6E6E',
'#6B6B6B',
'#696969',
'#666666',
'#636363',
'#616161',
'#5E5E5E',
'#5C5C5C',
'#595959',
'#575757',
'#545454',
'#525252',
'#4F4F4F',
'#4D4D4D',
'#4A4A4A',
'#474747',
'#454545',
'#424242',
'#404040',
'#3D3D3D',
'#3B3B3B',
'#383838',
'#363636',
'#333333',
'#303030',
'#2E2E2E',
'#2B2B2B',
'#292929',
'#262626',
'#242424',
'#212121',
'#1F1F1F',
'#1C1C1C',
'#1A1A1A',
'#171717',
'#141414',
'#121212',
'#0F0F0F',
'#0D0D0D',
'#0A0A0A',
'#080808',
'#050505',
'#030303']

# rainbowcolors=['#3A5FCD','#87CEFA', '#00CED1','#D02090','#DC143C']
rainbowcolors=['dodgerblue','darkseagreen',"orange", "salmon", '#D02090','#DC143C']


'''zbin=0.3
lowradius = radius[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20) &(mass>10.5)].flatten()
lowmass = mass[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20)&(mass>10.5)].flatten()
lowfit = np.polyfit(lowmass, lowradius, 2)
lowquad,lowlin,lowinter=np.poly1d(lowfit)
lowbestfit = lowinter + lowlin * (lowmass ) + lowquad * (lowmass) ** 2

zbin=0.45
middleradius = radius[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth))  & (radius > 1) & (radius<20)&(mass>10.5)].flatten()
middlemass = mass[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20)&(mass>10.5)].flatten()
middlefit = np.polyfit(middlemass, middleradius, 2)
midquad,midlin,midinter=np.poly1d(middlefit)
midbestfit = midinter + midlin * (middlemass ) + midquad * (middlemass) ** 2

zbin=0.65
highradius = radius[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth))  & (radius > 1) & (radius<20)&(mass>10.5)].flatten()
highmass = mass[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth))  & (radius > 1) & (radius<20)&(mass>10.5)].flatten()
highfit = np.polyfit(highmass, highradius, 2)
highquad,highlin,highinter=np.poly1d(highfit)
highbestfit = highinter + highlin * (highmass ) + highquad * (highmass) ** 2'''

# fig, ax = plt.subplots(1)
plt.figure(facecolor="white")
'''plt.scatter( lowradius,lowmass, label=("0.3"), linewidth=0.75, alpha=0.3, color="mediumpurple", marker="*", edgecolor=None)
plt.scatter(  middleradius,middlemass, label=("0.45"), , linewidth=0.75, alpha=0.2, color="green", marker="*", edgecolor=None)
plt.scatter( highradius,highmass, label=("0.65"), linewidth=0.75, alpha=0.1, color="orangered", marker="*", edgecolor=None)
'''
k=0
print np.arange(0,len(emceeredshifts),10)

for i in np.arange(0,len(emceeredshifts),10):
	zbin=emceeredshifts[i]
# 	print zbin
	lowradius = radius[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20) &(mass>10.5)].flatten()
# 	print lowradius
	lowmass = mass[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20)&(mass>10.5)].flatten()
	if len(lowradius)>0:
		lowfit = np.polyfit(lowmass, lowradius, 1)
		lowlin,lowinter=np.poly1d(lowfit)
		lowbestfit = lowinter + lowlin * (lowmass )
		ind = np.argsort(lowmass)
# 		plt.scatter( lowmass, lowbestfit, color=rainbowcolors[k], marker='.', label=str(zbin))
		plt.plot( lowmass[ind], lowbestfit[ind], color=rainbowcolors[k], label=str(zbin), linewidth=4)
# 	plt.scatter( lowmass, lowbestfit, color=rainbowcolors[k], marker='.', edgecolor='none', label=str(zbin))
	ind=np.argsort(lowmass)
	plt.scatter( lowmass[ind], lowradius[ind], color=rainbowcolors[k], marker='.', alpha=.2, edgecolor='none')
	k+=1
	# lowfit = np.polyfit(lowmass, lowradius, 2)
# 	lowquad,lowlin,lowinter=np.poly1d(lowfit)
# 	lowbestfit = lowinter + lowlin * (lowmass ) + lowquad * (lowmass) ** 2
# 	plt.scatter( lowmass, lowbestfit, color=grayscalecolors[i], marker='.', edgecolor='none', label=str(zbin))
	

# plt.plot( middlemass, midbestfit, color="lightpink", label=("0.45"))
# plt.plot( highmass, highbestfit, color="orangered", label=("0.65"))

plt.legend(loc='lower right')

plt.title("Galaxy radius vs mass\nquadratic mcmc mass predictions")
plt.xlabel("$Log$ $M_\odot$")
plt.ylabel("$Physical$ $radius$ $(kpc)$")
plt.ylim(2,15)
plt.xlim(10.6,11.8)

plt.yscale('log')
plt.show()
	