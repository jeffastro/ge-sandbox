import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from scipy import stats
from astropy.cosmology import FlatLambdaCDM
from scipy.optimize import minimize

def powerlawcomponents(xdata, ydata):
        lnx = np.log(xdata)
        lny = np.log(ydata)
        slope, intercept, r, p, std_err = stats.linregress(lnx, lny)
        bestfit = (np.exp(intercept))*xdata**(slope)
        return [slope, intercept]


def newpowerlawcomponents (xdata, ydata,pivot):
	lnx = np.log(xdata) - np.log(pivot)
	lny = np.log(ydata)
	slope, intercept, r_value, p_value, std_err = stats.linregress(lnx,lny)
	bestfit = (np.exp(intercept))*xdata**(slope)
	return [slope, intercept]

def powerlawfit (xdata, ydata):
	lnx = np.log(xdata)
	lny = np.log(ydata)
	slope, intercept, r_value, p_value, std_err = stats.linregress(lnx,lny)
	bestfit = (np.exp(intercept))*xdata**(slope)
	return bestfit

def powerline (params, xdata):
	a,b,c = params
	return (a+b*np.log(1+xdata))

def loglike(params, xdata, ydata):
	EPSILON = np.power(2., -16)
	sigma=params[2]
	#--------ln Likelihood--------
	yPred = ((np.log(ydata) - powerline(params, xdata))**2 / (2*((EPSILON+sigma)**2)))
	# Calculate negative log likelihood
# 	LL = -np.sum(stats.norm.logpdf(ydata, loc=yPred, scale=sigma))
	LL = np.sum(yPred)
	return(LL)


cosmo = FlatLambdaCDM(H0 =70, Om0 = 0.3)


thecat = fits.getdata('matched_cs82_redmag_040316.fits')
emceedata = fits.getdata("emceeRun051716_quad.fits")

intercepts = emceedata['a'].copy()
lins = emceedata['b'].copy()
quads = emceedata['c'].copy()

mass_int_err = emceedata['sigma'].copy()
factors = emceedata['factor'].copy()
emceeredshifts = emceedata['redshift'].copy()
galnumbers = emceedata['number_of_galaxies'].copy()


spher = thecat['SPHEROID_REFF_WORLD'].copy()
redshift = thecat['ZREDMAGIC'].copy()
zband = thecat['ZBAND'].copy()



#get pivot point
kavli_fits = fits.getdata("matchedkavli051716.fits")
zbandforpivot = kavli_fits['ZBAND'].copy()
pivot = np.median(zbandforpivot)
print pivot


radius = []
mass = np.zeros(len(zband)).tolist()


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

# -------fits  power laws through 3 mass bins and plots.------------
radii_11=[]
errors_11=[]
radii_11_2=[]
errors_11_2=[]
radii_11_4=[]
errors_11_4=[]

currentmstar = 11
for n in xrange(len(emceeredshifts)):
    zbin = emceeredshifts[n]
    currentinterr = mass_int_err[n]
    currentradius = radius[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20) ].flatten()
    massofinterest = mass[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20)].flatten()
    currentfit = np.polyfit(massofinterest, currentradius, 2)
    currentquad,currentlin,currentinter=np.poly1d(currentfit)
    bestfitradius = currentinter + currentlin * (currentmstar ) + currentquad * (currentmstar) ** 2
    radii_11.append(bestfitradius)
    bestfiterror =  currentlin * (currentinterr ) + currentquad * (currentinterr) ** 2
    print "current intrinsic scatter: " + str(currentinterr)
    print "current bestfit error: " + str(bestfiterror)
    print currentinter
    print currentlin
    print currentquad
    errors_11.append(bestfiterror)

currentmstar = 11.2
for n in xrange(len(emceeredshifts)):
    zbin = emceeredshifts[n]
    currentinterr = mass_int_err[n]
    currentradius = radius[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20)].flatten()
    massofinterest = mass[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20)].flatten()
    currentfit = np.polyfit(massofinterest, currentradius, 2)
    currentquad,currentlin,currentinter=np.poly1d(currentfit)
    bestfitradius = currentinter + currentlin * (currentmstar ) + currentquad * (currentmstar) ** 2
    radii_11_2.append(bestfitradius)
    bestfiterror = currentlin * (currentinterr ) + currentquad * (currentinterr) ** 2
    errors_11_2.append(bestfiterror)
	
currentmstar = 11.4
for n in xrange(len(emceeredshifts)):
    zbin = emceeredshifts[n]
    currentinterr = mass_int_err[n]
    currentradius = radius[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20) ].flatten()
    massofinterest = mass[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20)].flatten()
    currentfit = np.polyfit(massofinterest, currentradius, 2)
    currentquad,currentlin,currentinter=np.poly1d(currentfit)
    bestfitradius = currentinter + currentlin * (currentmstar ) + currentquad * (currentmstar) ** 2
    radii_11_4.append(bestfitradius)
    bestfiterror =  currentlin * (currentinterr ) + currentquad * (currentinterr) ** 2
    errors_11_4.append(bestfiterror)


# bestlowfit = powerlawfit(emceeredshifts, radii_11)
# bestmidfit = powerlawfit(emceeredshifts, radii_11_2)
# besthighfit = powerlawfit(emceeredshifts, radii_11_4)

radii_11 = np.array(radii_11)
initParams = powerlawcomponents(emceeredshifts, radii_11)
initParams.append(0.1)

results = minimize(loglike, initParams, args=(emceeredshifts, radii_11), method='BFGS')

estParms = results.x
 
intercept = estParms[0]
slope = estParms[1]
trylowfit = powerlawfit((1+emceeredshifts), radii_11)
trymedfit = powerlawfit((1+emceeredshifts), radii_11_2)
tryhighfit = powerlawfit((1+emceeredshifts), radii_11_4)
	
plt.figure(facecolor='white')

plt.plot(emceeredshifts, trylowfit, color="darkseagreen")
plt.plot(emceeredshifts, trymedfit, color="dodgerblue")
plt.plot(emceeredshifts, tryhighfit, color="darkorchid")

plt.errorbar(emceeredshifts, radii_11, yerr=errors_11, marker='.', mec='none',mfc='darkseagreen',linestyle="None", label="$11.0$ $Log$ $M_\odot$")
plt.errorbar(emceeredshifts, radii_11_2, yerr=errors_11_2, marker='.', mec='none',mfc='dodgerblue', label="$11.2$ $Log$ $M_\odot$", linestyle="None")
plt.errorbar(emceeredshifts, radii_11_4, yerr=errors_11_4, marker='.', mec='none',mfc='darkorchid', label="$11.4$ $Log$ $M_\odot$", linestyle="None")

plt.legend(loc="upper right", fontsize=8)
plt.ylabel("$Physical$ $radius$ $(kpc)$")
plt.xlabel("$Redshift$")
plt.title("$Radius$ $vs$ $redshift$ $in$ $stellar$ $mass$ $bins$\n$power$ $law$ $fit$ $using$ $(1+z)$")
plt.xlim([0.2,0.7])
plt.show()


currentmassbinradii = []
slopes = []
mstarbins = np.linspace(11.0, 11.5, 10)
for currentmstar in mstarbins:
        currentmassbinradii=[]
        for zbin in emceeredshifts:
                currentradius = radius[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20) ].flatten()
                massofinterest = mass[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (radius > 1) & (radius<20)].flatten()
                currentfit = np.polyfit(massofinterest, currentradius, 2)
                currentquad,currentlin,currentinter=np.poly1d(currentfit)
                bestfitradius = currentinter + currentlin * (currentmstar ) + currentquad * (currentmstar) ** 2
                currentmassbinradii.append(bestfitradius)
        currentmassbinradii = np.array(currentmassbinradii)
        print (emceeredshifts)
        print emceeredshifts.shape
        print(currentmassbinradii)
        print currentmassbinradii.shape
        print(mstarbins)
        bestslope, bestintercept = newpowerlawcomponents((1+emceeredshifts), currentmassbinradii, 0.45)
        slopes.append(bestslope)

plt.figure(facecolor='white')
plt.scatter(mstarbins,slopes, color = "darkorchid", edgecolor = 'none')
plt.title("$Power$ $law$ $exponential$ $coefficients$")
plt.ylabel("$exponential$ $coefficient$ $on$ $radius$ $vs$ $redshift$")
plt.xlabel("$M_\odot$")
plt.show()


'''
-----attempt at power law with (1+z)------



import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from scipy import stats
from astropy.cosmology import FlatLambdaCDM

def powerlawfit (xdata, ydata):
	lnx = np.log(xdata)
	lny = np.log(ydata)
	slope, intercept, r_value, p_value, std_err = stats.linregress(lnx,lny)
	bestfit = (np.exp(intercept))*xdata**(slope)
	return bestfit

def powerlawcomponents (xdata, ydata):
	lnx = np.log(xdata)
	lny = np.log(ydata)
	slope, intercept, r_value, p_value, std_err = stats.linregress(lnx,lny)
	bestfit = (np.exp(intercept))*xdata**(slope)
	return slope, intercept

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


radius = []
mass = np.zeros(len(zband)).tolist()


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
	
oneplusz = [1+z for z in emceeredshifts]
'''
# slope, intercept = powerlaw components(oneplusz, radii_10_8)
# bestlowfit = (np.exp(intercept))*emceeredshifts**(slope)
# slope, intercept = powerlawcomponents(oneplusz, radii_11_1)
# bestmidfit = (np.exp(intercept))*emceeredshifts**(slope)
# slope, intercept = powerlawcomponents(oneplusz, radii_11_4)
# besthighfit = (np.exp(intercept))*emceeredshifts**(slope)
'''

slope, intercept = powerlawcomponents(oneplusz, radii_10_8)
bestlowfit = (np.exp(intercept))*oneplusz**(slope)
slope, intercept = powerlawcomponents(oneplusz, radii_11_1)
bestmidfit = (np.exp(intercept))*oneplusz**(slope)
slope, intercept = powerlawcomponents(oneplusz, radii_11_4)
besthighfit = (np.exp(intercept))*oneplusz**(slope)


	
plt.figure(facecolor='white')
plt.scatter(emceeredshifts, radii_10_8, marker='o', edgecolor='none',color='darkseagreen', label="$10.8$ $Log$ $M_\odot$")
plt.scatter(emceeredshifts, radii_11_1, marker='o', edgecolor='none',color='dodgerblue', label="$11.1$ $Log$ $M_\odot$")
plt.scatter(emceeredshifts, radii_11_4, marker='o', edgecolor='none',color='darkorchid', label="$11.4$ $Log$ $M_\odot$")
# plt.plot(emceeredshifts, bestlowfit,color='darkseagreen')
# plt.plot(emceeredshifts, bestmidfit,color='dodgerblue')
# plt.plot(emceeredshifts, besthighfit,color='darkorchid')
plt.legend(loc="lower left")
plt.ylabel("$Physical$ $radius$ $(kpc)$")
plt.xlabel("$1+Redshift$")
plt.title("$Radius$ $vs$ $redshift$ $in$ $stellar$ $mass$ $bins$\n$power$ $law$ $fit$")
plt.xlim([0.2,0.7])
plt.show()



'''
