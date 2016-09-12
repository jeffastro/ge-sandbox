# Last modified 05.17.16 at 00:28
# loop for creating many linear fits across z bins with a single pivot for all.
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')


import numpy as np
np.seterr(invalid='ignore')
import matplotlib.pyplot as plt
from astropy.io import fits
import corner
import emcee
from datetime import datetime
import sys



# -------Ln Likelihood function----------------------------------------
def lnlike(params, x, y, emass, ezband, pivot):
    print
    print
    print
    print "DERIVE THE ERROR EQUATION FIRST"
    sys.exit()
    a, b, c, sigmaint, factor = params
    sigma = np.sqrt(sigmaint**2 + emass**2 + (c*ezband)**2)
    model = a + b * (x - pivot) + c * (x - pivot) ** 2
    inthesum = np.log(
        factor * (1 / ((np.sqrt(2 * np.pi)) * sigma)) * (np.exp(-(np.square(y - model)) / (2 * (sigma ** 2)))) \
        + 0.625 * (1 - factor))
    return np.sum(inthesum)


def linlnlike(params, x, y, emass, ezband, pivot):
    a, b, sigmaint, factor = params
    sigma = np.sqrt(np.abs(sigmaint**2 + emass**2 + b*(ezband)**2))
    model = a + b * (x - pivot) 
    inthesum = factor * (1 / ((np.sqrt(2 * np.pi)) * sigma)) * (np.exp(-(np.square(y - model)) / (2 * (sigma ** 2)))) + 0.625 * (1 - factor)
    if np.sum(inthesum)<0.0:
        return -np.inf
    return np.sum(np.log(inthesum))


def linlnprior(params):
    a, b, sigmaint, factor = params
    if 0 < sigmaint <= 1 and 0.5 < factor < 1 and 5 < a < 45 and -6 < b < 0:
        return 0.0
    return -np.inf


def lnprior(params):
    a, b, c, sigmaint, factor = params
    if 0 < sigmaint <= 1 and 0.5 < factor < 1 and 5 < a < 45 and -6 < b < 0 and -1 < c < 2:
        return 0.0
    return -np.inf


def lnprob(params, x, y, emass, ezband, pivot, linearmode=True):
    if linearmode:
        lp = linlnprior(params)
        if not np.isfinite(lp):
            return -np.inf
        return lp + linlnlike(params, x, y, emass, ezband, pivot)
    else:
        print
        print
        print
        print "DERIVE THE ERROR EQUATION FIRST"
        sys.exit()
        lp = lnprior(params)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike(params, x, y, emass, ezband, pivot)

def todaystring():
    now = datetime.now()
    return "0" + str(now.month)+str(now.day)+str(now.year)

# ------get data-------------------------------------------------------

# catpath = raw_input("catalog path and filename? ")
# mode = raw_input("linear or quadratic? ").lower()

catpath = "matchedkavli051716.fits"
mode = "lin"

if (mode == "linear") or (mode == "lin") or (mode == "l"): linearmode = True
else: linearmode = False

kavli_fits = fits.getdata(catpath)
mstar = kavli_fits['MASS_IR_BEST'].copy()
zband = kavli_fits['ZBAND'].copy()
redshift = kavli_fits['ZREDMAGIC'].copy()
emass = kavli_fits['MASSERR_IR_BEST'].copy()
ezband = kavli_fits['ZBANDERR'].copy()

for i in xrange(len(emass)):
    if emass[i]>0.0:
        pass
    else:
        emass[i]=0.0
        
for i in xrange(len(ezband)):
    if ezband[i]>0.0:
        pass
    else:
        ezband[i]=0.0

n = -1
'''n = input("Speaking of redshift, we can start at any redshift in increments of \
0.02. Entering 0 begins the chain at z=0.24, 1 puts you at z=0.26, and so on. What \
number are we working with today? ")'''
zbin = 0.24 + n*0.02
width = 0.02
halfwidth = width / 2
pivot = np.median(zband)

a = []
alow = []
ahigh = []
b = []
blow = []
bhigh = []
if linearmode:
    c = []
    clow = []
    chigh = []
sigmas = []
siglow = []
sighigh = []
factors = []
factorlow = []
factorhigh = []
redshifts = []
galnumbers = []

# csv=True
csv = False


if n>-1:
    limitontheoperation = zbin + halfwidth
    csv = True
else:
    limitontheoperation = np.max(redshift)
    csv = False
    n=1

while (zbin < limitontheoperation):
    zbandofinterest = zband[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin -
                             halfwidth)) & (mstar > 0) & (zband > 0)]
    lowestz = np.percentile(zbandofinterest, 99)
    medz_mstar = mstar[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth))
                        & (mstar > 0) & (zband > 0) & (zband < lowestz)]
    medz_zband = zband[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (mstar > 0) & (zband > 0) & (zband < lowestz)]
    emasscurrent = emass[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (mstar > 0) & (zband > 0) & (zband < lowestz)]
    ezbandcurrent = ezband[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (mstar > 0) & (zband > 0) & (zband < lowestz)]
    
    ydata = medz_mstar
    xdata = medz_zband
    
    # -------Markov chain---------------------------------------------------
    nwalkers = 10
    ndim = 4
    if linearmode==False:
        ndim+=1
    if linearmode:
        pos = np.array([10 + (np.random.ranf(nwalkers) * 10), -3 + np.random.ranf(nwalkers), np.random.ranf(nwalkers), 0.5 + (np.random.ranf(nwalkers) / 2.1)])
    else:
        pos = np.array([36 + (np.random.ranf(nwalkers) * 10), -3 + np.random.ranf(nwalkers), (np.random.ranf(nwalkers) / 10),np.random.ranf(nwalkers), 0.5 + (np.random.ranf(nwalkers) / 2.1)])
    pos = pos.transpose()
    pos = pos.tolist()
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(xdata, ydata, emasscurrent, ezbandcurrent, pivot, linearmode))
    ntrials = 10000
    sampler.run_mcmc(pos, ntrials)
    burnin = 5000
    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
    lnprobs = sampler.flatlnprobability
    if linearmode:
        samples = samples[samples[:, 3] > 0.6]
    else:
        samples = samples[samples[:, 4] > 0.6]
    
    # -------Save the sample results-----------------------------------------
    adata = sampler.flatchain[burnin:, 0]
    bdata = sampler.flatchain[burnin:, 1]
    if linearmode:
        sigmadata = sampler.flatchain[burnin:, 2]
        factordata = sampler.flatchain[burnin:, 3]
    else:
        cdata = sampler.flatchain[burnin:, 2]
        sigmadata = sampler.flatchain[burnin:, 3]
        factordata = sampler.flatchain[burnin:, 4]
    adata = adata[factordata > 0.6]
    bdata = bdata[factordata > 0.6]
    if linearmode == False:
        cdata = cdata[factordata > 0.6]
    sigmadata = sigmadata[factordata > 0.6]
    factordata = factordata[factordata > 0.6]
    
    # -------Corner plot-----------------------------------------------------
    if linearmode:
        figure = corner.corner(samples, labels=["$a$", "$b$", "$\sigma_{int}$", "$factor$"], quantiles=[0.16, 0.5, 0.84], show_titles=True)
    else:
        figure = corner.corner(samples, labels=["$a$", "$b$", "$c$", "$\sigma_{int}$", "$factor$"], quantiles=[0.16, 0.5, 0.84], show_titles=True)
    figure.gca().annotate(str(mode) + "MSTAR vs ZBAND model,\nredshift = " + str(zbin) +
                          u"\u00B1" + str(halfwidth), xy=(0.5, 1.0),
                          xycoords="figure fraction", xytext=(0, -5),
                          textcoords="offset points", ha="center", va="top")
    figure.savefig("MarkovChainPDF"+todaystring()+"_" + str(mode)+"_{0}.png".format(n))

    # -------Compute quantiles----------------------------------------
    a_low = np.percentile(adata, [16])
    a_mid = np.percentile(adata, [50])
    a_high = np.percentile(adata, [84])
    a_low = a_low[0]
    a_mid = a_mid[0]
    a_high = a_high[0]
    a_mcmc = [a_mid, a_high - a_mid, a_mid - a_low]

    b_low = np.percentile(bdata, [16])
    b_mid = np.percentile(bdata, [50])
    b_high = np.percentile(bdata, [84])
    b_low = b_low[0]
    b_mid = b_mid[0]
    b_high = b_high[0]
    b_mcmc = [b_mid, b_high - b_mid, b_mid - b_low]
    
    if linearmode==False:
        c_low = np.percentile(cdata, [16])
        c_mid = np.percentile(cdata, [50])
        c_high = np.percentile(cdata, [84])
        c_low = c_low[0]
        c_mid = c_mid[0]
        c_high = c_high[0]
        c_mcmc = [c_mid, c_high - c_mid, c_mid - c_low]

    s_low = np.percentile(sigmadata, [16])
    s_mid = np.percentile(sigmadata, [50])
    s_high = np.percentile(sigmadata, [84])
    s_low = s_low[0]
    s_mid = s_mid[0]
    s_high = s_high[0]
    s_mcmc = [s_mid, s_high - s_mid, s_mid - s_low]

    f_low = np.percentile(factordata, [16])
    f_mid = np.percentile(factordata, [50])
    f_high = np.percentile(factordata, [84])
    f_low = f_low[0]
    f_mid = f_mid[0]
    f_high = f_high[0]
    f_mcmc = [f_mid, f_high - f_mid, f_mid - f_low]

    # ------Add info to arrays for the final fits file------------
    a.append(a_mid)
    alow.append(a_mcmc[1])
    ahigh.append(a_mcmc[2])
    b.append(b_mid)
    blow.append(b_mcmc[1])
    bhigh.append(b_mcmc[2])
    if linearmode==False:
        c.append(c_mid)
        clow.append(c_mcmc[1])
        chigh.append(c_mcmc[2])
    sigmas.append(s_mid)
    siglow.append(s_mcmc[1])
    sighigh.append(s_mcmc[2])
    factors.append(f_mid)
    factorlow.append(f_mcmc[1])
    factorhigh.append(f_mcmc[2])
    redshifts.append(zbin)
    galnumbers.append(len(medz_mstar))

    if linearmode:
        annotation = """MCMC result:
$a = {0[0]:5.3f}+{0[1]:5.3f}-{0[2]:5.3f}$
$b = {1[0]:5.3f}+{1[1]:5.3f}-{1[2]:5.3f}$
$\sigma_int$ $=\ {2[0]:5.3f}+{2[1]:5.3f}-{2[2]:5.3f}$
$factor = {5[0]:5.3f}+{5[1]:5.3f}-{5[2]:5.3f}$
${3}$ $walkers,$ ${4}$ $steps,$ ${6}$ $discarded$""".format(a_mcmc, b_mcmc, s_mcmc,
                                                            nwalkers, ntrials, f_mcmc, burnin)
    else:
        annotation = """MCMC result:
$a = {0[0]:5.3f}+{0[1]:5.3f}-{0[2]:5.3f}$
$b = {1[0]:5.3f}+{1[1]:5.3f}-{1[2]:5.3f}$
$c = {7[0]:5.3f}+{7[1]:5.3f}-{7[2]:5.3f}$
$\sigma_i$ = ${2[0]:5.3f}+{2[1]:5.3f}-{2[2]:5.3f}$
$factor = {5[0]:5.3f}+{5[1]:5.3f}-{5[2]:5.3f}$
${3}$ $walkers,$ ${4}$ $steps,$ ${6}$ $discarded$""".format(a_mcmc, b_mcmc, s_mcmc,
                                                            nwalkers, ntrials, f_mcmc,
                                                            burnin, c_mcmc)
    plt.figure()
    xranges = np.arange(np.min(medz_zband) - 0.2, np.max(medz_zband) + 0.2, 0.1)

    if linearmode:
        bestfit = a_mcmc[0] + b_mcmc[0] * (xranges -pivot)
    else:
        bestfit = a_mcmc[0] + b_mcmc[0] * (xranges - pivot) + c_mcmc[0] * \
                  (xranges - pivot) ** 2

    plt.plot(xranges, bestfit, color="indigo")
    plt.errorbar(medz_zband, medz_mstar, yerr=emasscurrent, xerr=ezbandcurrent,
                 color='palevioletred', marker=".",capsize=0, linestyle="None",
                 alpha=0.2, label=annotation)
    plt.title('Galaxy M-star (MGC) v s ZBAND (stripe82)\nredshift = ' + str(zbin)
              + u"\u00B1" + str(halfwidth))
    plt.legend(loc='best', fontsize=8)
    plt.xlabel("ZBAND")
    plt.ylabel("$M*\ (Log\ M_\odot)$")
    plt.savefig("emcee_mstarzband_"+str(mode)+"_" + todaystring() + "_" + str(n) + ".png")
    plt.show()

    plt.figure()
    steps = np.arange(len(lnprobs))
    plt.plot(steps, lnprobs, color = "slateblue")
    plt.title("Log likehood in redshift = " + str(zbin) + u"\u00B1" + str(halfwidth))
    plt.xlabel("step")
    plt.ylabel("Ln $\mathcal{L}$")
    plt.ylim(np.max(lnprobs) - 30, np.max(lnprobs) + 5)
    plt.savefig("lnprobs" + todaystring() + "_" + str(n) + ".png")
    plt.show()

    zbin += 0.02
    n += 1
    print "yet another trial done"

results=[np.array(a), np.array(b), np.array(sigmas), np.array(factors),
         np.array(redshifts), np.array(galnumbers), np.array(siglow), np.array(sighigh),
         np.array(alow), np.array(ahigh), np.array(blow), np.array(bhigh),
         np.array(factorlow), np.array(factorhigh)]

if linearmode == False:
    results.append(np.array(c))
    results.append(np.array(clow))
    results.append(np.array(chigh))

finalheaders = ['a', 'b', 'sigma', 'factor', 'redshift', 'number_of_galaxies', 'sigmalow',
                'sigmahigh','a_low', 'a_high', 'b_low', 'b_high', 'factor_low',
                'factor_high']

columns = []

if linearmode == False:
    finalheaders.append('c') 
    finalheaders.append('c_low')
    finalheaders.append('c_high')

if csv:
    target = open('emceelin052516.csv', 'w')
    results = [a,b,sigmas,factors,redshift,galnumbers,siglow,sighigh,alow,\
ahigh,blow, bhigh, factorlow, factorhigh]
    for i in results:
        target.write(str(i[0]))
        target.write(',')
    target.write("\n")
    target.close()
else:
    for i in xrange(len(finalheaders)):
        columns.append(fits.Column(name=finalheaders[i],
                                    format='E', array=results[i]))
    
    tbhdu = fits.BinTableHDU.from_columns(columns)
    tbhdu.writeto('emceeRun' + todaystring() + '_' + str(mode) + '.fits')