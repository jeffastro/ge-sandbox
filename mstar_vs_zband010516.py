#Last edited 01.20.16 at 20:21
#Just doing one redshift bin and adding a fit line

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
from scipy import stats
from scipy.optimize import minimize
import pylab as py


kavli_fits = fits.getdata("matchedkavli112315.fits")
mstar = kavli_fits['MASS_IR_BEST'].copy()
zband = kavli_fits['ZBAND'].copy()
redshift = kavli_fits['ZREDMAGIC'].copy()
mstarerr = kavli_fits['MASSERR_IR_BEST'].copy()

medz_mstar=mstar[(redshift<=(0.501)) & (redshift>=(0.499)) & (mstar>0)& (zband>0)& (mstarerr>0)]                      
medz_zband=zband[(redshift<=(0.501)) & (redshift>=(0.499)) & (mstar>0)& (zband>0)& (mstarerr>0)]
medz_err=mstarerr[(redshift<=(0.501)) & (redshift>=(0.499)) & (mstar>0)& (zband>0)& (mstarerr>0)]


#-----------------Max Likelihood function----------

ydata = medz_mstar
xdata = medz_zband

def gaussian(params):
    k = params[0]
    b = params[1]   
    sigma = params[2]


#--------ln Likelihood--------
    yPred = (np.log(1 / (((2*np.pi)**0.5)*sigma))) - ((ydata - k*xdata - b)**2 / (2*(sigma**2)))

    # Calculate negative log likelihood
    LL = -np.sum( stats.norm.logpdf(ydata, loc=yPred, scale=sigma ) )

    return(LL)


#-----------------Fit line-------------------------


#init_vals = [1, 0, 1]
#best_vals, covar = curve_fit(gaussian, medz_zband, medz_mstar, sigma=medz_err)
#, sigma = medz_err, p0=(1./np.std(medz_mstar), np.argmax(medz_mstar) ,0,0,1))
#try:
#    best_vals, covar = curve_fit(gaussian, medz_zband, medz_mstar, p0 = None, sigma = None)

#try:
#    best_vals, covar = curve_fit(linear, medz_zband, medz_mstar)
#
#except RuntimeError:
#    print("Error - curve_fit failed")
    
initParams = [-.177, 14.525, 1]

results = minimize(gaussian, initParams, method='Nelder-Mead')
print results.x

estParms = results.x
 
k = estParms[0]
b = estParms[1]
sigma = estParms[2]
      
yOut = estParms[0] * xdata + estParms[1]



plt.scatter(medz_zband, medz_mstar, color='c', marker="*", label = 'z~0.5')
plt.plot(xdata, yOut, color = 'purple')#,best_vals[2]))
#plt.axis([15, 23, 9.5,12.5])
plt.title('Galaxy M-star (MGC) v s ZBAND (stripe82)')
plt.annotate('MASS = 14.481 - 0.174 * ZBAND\nIntrinsic scatter = 11.002', xy=(21,10.83), xytext=(20.55,11.2),\
arrowprops=dict(facecolor='black',width=2, shrink=0.05))
plt.xlabel("ZBAND")
plt.ylabel("MASS_IR_BEST")
##plt.xscale('log', basex=10)
##plt.yscale('log', basey=10)
#plt.legend(loc="best")
plt.show()
# plt.savefig("mstar_vs_zband_fit012016.pdf")