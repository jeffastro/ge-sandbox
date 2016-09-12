# -*- coding: utf-8 -*-
#Sasha Safonova


#Uses a Bayesian likelihood formulation with Markov chain monte carlo (mcmc)
#to estimate linear model parameters.
#Repeats the mcmc in redshift bins of 0.02 across an entire catalog, provided
#each bin has at least 100 galaxies.


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import corner #get it here https://github.com/dfm/corner.py
import emcee #get it here https://github.com/dfm/emcee


#-------Ln Likelihood function----------------------------------------
def lnlike(params, x, y):
    
    #k = slope
    #b = intercept
    #sigma = intrinsic scatter
    #factor = portion of data that is good
    k, b, sigma, factor = params
    model = k * x + b
    return np.sum(np.log(factor*(1 / ((np.sqrt(2*np.pi))*sigma))*(np.exp(-(np.square(y - model)) / (2*(sigma**2))))\
    + 0.625*(1-factor)))
    
def lnprior(params):
    k, b, sigma, factor = params
    
    #these are the priors that were appropriate for my data. 
    #choose priors that make sense with your catalog.
    if 0<sigma<=1 and 0.5<factor<1 and -1<k<0 and 5<b<25:
        return 0.0
    return -np.inf

def lnprob(params, x, y):
    lp = lnprior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(params, x, y)
    
def markovchain(xdata, ydata, nwalkers, ntrials, burnin):
    ndim = 4
    pos = [[-0.3+(np.random.ranf(nwalkers)/10), 15.1+np.random.ranf(nwalkers), 0.11+(np.random.ranf(nwalkers)/10), 0.5+(np.random.ranf(nwalkers)/2.1)]]
    pos = np.array(pos)
    pos=pos[0]
    pos = pos.transpose()
    pos=pos.tolist()
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = (xdata, ydata))
    sampler.run_mcmc(pos, ntrials)
    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
    return samples
    

#------get data-------------------------------------------------------
filename = input("Enter the name of your catalog: ")
catalog = fits.getdata(filename)
mstar = catalog['MASS_IR_BEST'].copy()
zband = catalog['ZBAND'].copy()
redshift = catalog['ZREDMAGIC'].copy()
                                                 
#zbin=0.5
width = 0.02
halfwidth = width/2

zbin = np.min(redshift)
n=1

invalidbins = []
needmorework = False
slopes=[]
slopelow=[]
slopehigh=[]
intercepts=[]
interlow=[]
interhigh=[]
sigmas=[]
siglow=[]
sighigh=[]
factors=[]
factorlow=[]
factorhigh=[]
redshifts=[]
galnumbers=[]
         
while (zbin<np.max(redshift)):
    medz_mstar=mstar[(redshift<=(zbin+halfwidth)) & (redshift>=(zbin-halfwidth)) & (mstar>0)& (zband>0) ]                      
    medz_zband=zband[(redshift<=(zbin+halfwidth)) & (redshift>=(zbin-halfwidth)) & (mstar>0)& (zband>0) ]
    
    #checks if a bin has at least 100 galaxies
    if(len(medz_mstar)<100):
        invalidbins.append([zbin, len(medz_mstar)])
        needmorework = True
    
    else:
        ydata = medz_mstar
        xdata = medz_zband

        #-------Markov chain---------------------------------------------------   
        nwalkers = 30
        ntrials = 50000
        burnin = 2000      
        samples = markovchain(xdata, ydata, nwalkers, ntrials, burnin)

 
        #-------Corner plot-----------------------------------------------------
        figure = corner.corner(samples, labels=["$slope$", "$intercept$", "$\sigma$","$factor$"],quantiles=[0.16, 0.5, 0.84], show_titles=True)
        figure.gca().annotate("Linear model,\nredshift = " + str(zbin) + u"\u00B1"+str(halfwidth), xy=(0.5, 1.0), xycoords="figure fraction",
                            xytext=(0, -5), textcoords="offset points",
                            ha="center", va="top")
        figure.savefig("MarkovChain_"+str(n)+".pdf")
        
        #-------Compute quantiles----------------------------------------
        k_low = np.percentile(samples[0], [16])
        k_mid = np.percentile(samples[0], [50])
        k_high = np.percentile(samples[0], [84])
        k_mcmc = [k_mid[0], k_high[0]-k_mid[0], k_mid[0]-k_low[0]]
        
        b_low = np.percentile(samples[1], [16])
        b_mid = np.percentile(samples[1], [50])
        b_high = np.percentile(samples[1], [84])
        b_mcmc = [b_mid[0], b_high[0]-b_mid[0], b_mid[0]-b_low[0]]
        
        
        s_low = np.percentile(samples[2], [16])
        s_mid = np.percentile(samples[2], [50])
        s_high = np.percentile(samples[2], [84])
        s_mcmc = [s_mid[0], s_high[0]-s_mid[0], s_mid[0]-s_low[0]]
        
        f_low = np.percentile(samples[3], [16])
        f_mid = np.percentile(samples[3], [50])
        f_high = np.percentile(samples[3], [84])
        f_mcmc = [f_mid[0], f_high[0]-f_mid[0], f_mid[0]-f_low[0]]
        
        
        annotation = """MCMC result:
            $slope = {0[0]:5.3f}+{0[1]:5.3f}-{0[2]:5.3f}$
            $intercept = {1[0]:5.3f}+{1[1]:5.3f}-{1[2]:5.3f}$
            $\sigma$ = ${2[0]:5.3f}+{2[1]:5.3f}-{2[2]:5.3f}$
            $factor = {5[0]:5.3f}+{5[1]:5.3f}-{5[2]:5.3f}$
            ${3}$ $walkers,$ ${4}$ $steps,$ ${6}$ $discarded$
        """.format(k_mcmc, b_mcmc,s_mcmc, nwalkers, ntrials, f_mcmc, burnin)
        
        slopes.append(k_mid)
        slopelow.append(k_mcmc[1])
        slopehigh.append(k_mcmc[2])
        intercepts.append(b_mid)
        interlow.append(b_mcmc[1])
        interhigh.append(b_mcmc[2])
        sigmas.append(s_mid)
        siglow.append(s_mcmc[1])
        sighigh.append(s_mcmc[2])
        factors.append(f_mid)
        factorlow.append(f_mcmc[1])
        factorhigh.append(f_mcmc[2])
        redshifts.append(zbin)
        galnumbers.append(len(medz_mstar))
        
        #-------Plot best fit and data---------------------------------------------
        plt.figure()
        xranges = np.arange(np.min(medz_zband)-0.2, np.max(medz_zband)+0.2, 0.1);
        bestfit = k_mcmc[0] * xranges + b_mcmc[0]
        plt.plot(xranges, bestfit, color="indigo")
        
        plt.scatter(medz_zband, medz_mstar, color='palevioletred', marker="*", alpha=0.8, label = annotation)
        plt.title('Galaxy M-star (MGC) v s ZBAND (stripe82)\nredshift = ' + str(zbin) + u"\u00B1"+str(halfwidth))
        plt.legend(loc='best', fontsize=8)
        plt.xlabel("ZBAND")
        plt.ylabel("MASS_IR_BEST")
        plt.show()
        plt.savefig("emcee_outliers_"+str(n)+".pdf")
    
    zbin=zbin+0.02
    n=n+1
if needmorework:
    np.savetxt("invalidbins.txt",invalidbins) 
    

c1=fits.Column(name = 'slope', format = 'E', array = np.array(slopes))
c2=fits.Column(name = 'intercept', format = 'E', array = np.array(intercepts))
c3=fits.Column(name = 'sigma', format = 'E', array = np.array(sigmas))
c4=fits.Column(name = 'factor', format = 'E', array = np.array(factors))
c5=fits.Column(name = 'redshift', format = 'E', array = np.array(redshifts))
c6=fits.Column(name = 'number_of_galaxies', format = 'E', array = np.array(galnumbers))
c7=fits.Column(name = 'sigma_low', format = 'E', array = np.array(siglow))
c8=fits.Column(name = 'sigma_high', format = 'E', array = np.array(sighigh))
c9=fits.Column(name = 'slope_low', format = 'E', array = np.array(slopelow))
c10=fits.Column(name = 'slope_high', format = 'E', array = np.array(slopehigh))
c11=fits.Column(name = 'intercept_low', format = 'E', array = np.array(interlow))
c12=fits.Column(name = 'intercept_high', format = 'E', array = np.array(interhigh))
c13=fits.Column(name = 'factor_low', format = 'E', array = np.array(factorlow))
c14=fits.Column(name = 'factor_high', format = 'E', array = np.array(factorhigh))

tbhdu = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14])
tbhdu.writeto('emceeRun.fits')
