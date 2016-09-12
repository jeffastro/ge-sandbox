# Last modified 05.17.16 at 00:28
# loop for creating many linear fits across z bins with a single pivot for all.
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import corner
import emcee


# import matplotlib.font_manager as font_manager

# -------Ln Likelihood function----------------------------------------
def lnlike(params, x, y, emass, ezband, pivot):
    a, b, c, sigmaint, factor = params
    sigma = np.sqrt(sigmaint**2 + emass**2 + (c*ezband)**2)
    model = a + b * (x - pivot) + c * (x - pivot) ** 2
    inthesum = np.log(
        factor * (1 / ((np.sqrt(2 * np.pi)) * sigma)) * (np.exp(-(np.square(y - model)) / (2 * (sigma ** 2)))) \
        + 0.625 * (1 - factor))

    return np.sum(inthesum)


def lnprior(params):
    a, b, c, sigmaint, factor = params
    if 0 < sigmaint <= 1 and 0.5 < factor < 1 and 5 < a < 45 and -6 < b < 0 and -1 < c < 2:
        return 0.0

    return -np.inf


def lnprob(params, x, y, emass, ezband, pivot):
    lp = lnprior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(params, x, y, emass, ezband, pivot)


# font_prop = font_manager.FontProperties(family="serif", size="x-small")

# ------get data-------------------------------------------------------
kavli_fits = fits.getdata("matchedkavli051716.fits")
mstar = kavli_fits['MASS_IR_BEST'].copy()
zband = kavli_fits['ZBAND'].copy()
redshift = kavli_fits['ZREDMAGIC'].copy()
emass = kavli_fits['MASSERR_IR_BEST'].copy()
ezband = kavli_fits['ZBANDERR'].copy()

n = 1
# zbin = np.min(redshift)
#zbin = 0.24
zbin = 0.24 + n*0.02
width = 0.02
halfwidth = width / 2
pivot = np.median(zband)

#n = 0
'''cute=range(1,41)
indices=[]
for item in cute:
    indices.append(str(item))'''

indices = ['1',
           '2',
           '3',
           '4',
           '5',
           '6',
           '7',
           '8',
           '9',
           '10',
           '11',
           '12',
           '13',
           '14',
           '15',
           '16',
           '17',
           '18',
           '19',
           '20',
           '21',
           '22',
           '23',
           '24',
           '25',
           '26',
           '27',
           '28',
           '29',
           '30',
           '31',
           '32',
           '33',
           '34',
           '35',
           '36',
           '37',
           '38',
           '39',
           '40']

a = []
alow = []
ahigh = []
b = []
blow = []
bhigh = []
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

# while (zbin<zbin+0.02):
while (zbin < np.max(redshift)):
    medz_mstar = mstar[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (mstar > 0) & (zband > 0)]
    medz_zband = zband[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (mstar > 0) & (zband > 0)]
    emasscurrent = emass[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (mstar > 0) & (zband > 0)]
    ezbandcurrent = ezband[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (mstar > 0) & (zband > 0)]
    
    ydata = medz_mstar
    xdata = medz_zband
    
    # -------Markov chain---------------------------------------------------
    nwalkers = 20
    ndim = 5
    pos = np.array([[36 + (np.random.ranf(nwalkers) * 10), -3 + np.random.ranf(nwalkers), (np.random.ranf(nwalkers) / 10),np.random.ranf(nwalkers), 0.5 + (np.random.ranf(nwalkers) / 2.1)]])
    
    #    pos = np.array(pos)
    pos = pos[0]
    pos = pos.transpose()
    pos = pos.tolist()
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(xdata, ydata, emasscurrent, ezbandcurrent, pivot))
    ntrials = 5000
    sampler.run_mcmc(pos, ntrials)
    burnin = 2000
    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
    samples = samples[samples[:, 4] > 0.6]
    
    # -------Save the sample results-----------------------------------------
    adata = sampler.flatchain[burnin:, 0]
    bdata = sampler.flatchain[burnin:, 1]
    cdata = sampler.flatchain[burnin:, 2]
    sigmadata = sampler.flatchain[burnin:, 3]
    factordata = sampler.flatchain[burnin:, 4]
    adata = adata[factordata > 0.6]
    bdata = bdata[factordata > 0.6]
    cdata = cdata[factordata > 0.6]
    sigmadata = sigmadata[factordata > 0.6]
    factordata = factordata[factordata > 0.6]
    
    # -------Corner plot-----------------------------------------------------
    figure = corner.corner(samples, labels=["$a$", "$b$", "$c$", "$\sigma_int$", "$factor$"], quantiles=[0.16, 0.5, 0.84], show_titles=True)
    figure.gca().annotate("Linear MSTAR vs ZBAND model,\nredshift = " + str(zbin) + u"\u00B1" + str(halfwidth),xy=(0.5, 1.0), xycoords="figure fraction",xytext=(0, -5), textcoords="offset points",ha="center", va="top")
    currentno = indices[n]
    figure.savefig("MarkovChainPDF051716_quad_{0}.pdf".format(currentno))

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

    annotation = """MCMC result:
        $a = {0[0]:5.3f}+{0[1]:5.3f}-{0[2]:5.3f}$
        $b = {1[0]:5.3f}+{1[1]:5.3f}-{1[2]:5.3f}$
        $c = {7[0]:5.3f}+{7[1]:5.3f}-{7[2]:5.3f}$
        $\sigma$ = ${2[0]:5.3f}+{2[1]:5.3f}-{2[2]:5.3f}$
        $factor = {5[0]:5.3f}+{5[1]:5.3f}-{5[2]:5.3f}$
        ${3}$ $walkers,$ ${4}$ $steps,$ ${6}$ $discarded$
    """.format(a_mcmc, b_mcmc, s_mcmc, nwalkers, ntrials, f_mcmc, burnin, c_mcmc)

    plt.figure()
    xranges = np.arange(np.min(medz_zband) - 0.2, np.max(medz_zband) + 0.2, 0.1);
    bestfit = a_mcmc[0] + b_mcmc[0] * (xranges - pivot) + c_mcmc[0] * (xranges - pivot) ** 2
    # bestfit = b_mcmc[0] + a_mcmc[0] * (xranges -pivot)

    plt.plot(xranges, bestfit, color="indigo")

    plt.scatter(medz_zband, medz_mstar, color='palevioletred', marker="*", alpha=0.8, label=annotation)
    plt.title('Galaxy M-star (MGC) v s ZBAND (stripe82)\nredshift = ' + str(zbin) + u"\u00B1" + str(halfwidth))
    plt.legend(loc='best', fontsize=8)
    plt.xlabel("ZBAND")
    plt.ylabel("MASS_IR_BEST")
    plt.show()
    plt.savefig("emcee_mstarzband_quad_051716_" + currentno + ".pdf")
    zbin += 0.02
    n += 1
    print "yet another trial done"

c1 = fits.Column(name='a', format='E', array=np.array(a))
c2 = fits.Column(name='b', format='E', array=np.array(b))
c3 = fits.Column(name='sigma', format='E', array=np.array(sigmas))
c4 = fits.Column(name='factor', format='E', array=np.array(factors))
c5 = fits.Column(name='redshift', format='E', array=np.array(redshifts))
c6 = fits.Column(name='number_of_galaxies', format='E', array=np.array(galnumbers))
c7 = fits.Column(name='sigmalow', format='E', array=np.array(siglow))
c8 = fits.Column(name='sigmahigh', format='E', array=np.array(sighigh))
c9 = fits.Column(name='a_low', format='E', array=np.array(alow))
c10 = fits.Column(name='a_high', format='E', array=np.array(ahigh))
c11 = fits.Column(name='b_low', format='E', array=np.array(blow))
c12 = fits.Column(name='b_high', format='E', array=np.array(bhigh))
c13 = fits.Column(name='factor_low', format='E', array=np.array(factorlow))
c14 = fits.Column(name='factor_high', format='E', array=np.array(factorhigh))
c15 = fits.Column(name='c', format='E', array=np.array(c))
c16 = fits.Column(name='c_low', format='E', array=np.array(clow))
c17 = fits.Column(name='c_high', format='E', array=np.array(chigh))


tbhdu = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17])

now = datetime.datetime.now()



tbhdu.writeto('emceeRun'+"0"+str(now.month)+str(now.day)+str(now.year)+'_quad.fits')
