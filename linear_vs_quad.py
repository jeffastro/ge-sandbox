from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
from scipy.stats import linregress

filename1="emceeRun031616_quad.fits"
quadata=fits.getdata(filename1)

filename2="emceeRun022816_linear.fits"
lindata=fits.getdata(filename2)


qa0= quadata['a'].copy()
qb0 = quadata['b'].copy()
qc0 = quadata['c'].copy()
qsigmas0 = quadata['sigma'].copy()
qfactors0 = quadata['factor'].copy()
qredshifts = quadata['redshift'].copy()
qgalnumbers = quadata['number_of_galaxies'].copy()
qsigmas1 = quadata['sigmalow'].copy()
qsigmas2 = quadata['sigmahigh'].copy()
qa1= quadata['a_low'].copy()
qb1 = quadata['b_low'].copy()
qc1 = quadata['c_low'].copy()
qa2= quadata['a_high'].copy()
qb2 = quadata['b_high'].copy()
qc2 = quadata['c_high'].copy()
qfactors1 = quadata['factor_low'].copy()
qfactors2 = quadata['factor_high'].copy()
qredshifts = quadata['redshift'].copy()


la0= lindata['a'].copy()
lb0 = lindata['b'].copy()
lsigmas0 = lindata['sigma'].copy()
lfactors0 = lindata['factor'].copy()
lredshifts = lindata['redshift'].copy()
lgalnumbers = lindata['number_of_galaxies'].copy()
lsigmas1 = lindata['sigmalow'].copy()
lsigmas2 = lindata['sigmahigh'].copy()
la1= lindata['a_low'].copy()
lb1 = lindata['b_low'].copy()
la2= lindata['a_high'].copy()
lb2 = lindata['b_high'].copy()
lfactors1 = lindata['factor_low'].copy()
lfactors2 = lindata['factor_high'].copy()
lredshifts = lindata['redshift'].copy()

plt.figure()
plt.suptitle("Comparison of slope and intercept in linear and quadratic parameter fits")
# plt.scatter(qredshifts, qfactors0, color="pink", label="quad")
# plt.scatter(lredshifts, lfactors0, color="blue", label="lin")


exes=np.linspace(np.min(lredshifts), np.max(lredshifts))
ax1 = plt.subplot2grid((2,1), (0,0))
ax2 = plt.subplot2grid((2,1), (1,0))

ax1.scatter(qredshifts, qa0, marker="s",color="skyblue", label = "quad")
ax1.scatter(lredshifts, lb0, marker="o",color="darkviolet",alpha =0.8,label = "lin")
ax1.set_ylabel("intercept")
#ax1.set_xlabel
ax1.legend(loc=4)
ax1.set_xlim([0.2,0.7])
plt.setp(ax1.get_xticklabels(), visible=False)

ax2.scatter(qredshifts, qb0, marker="s",color="skyblue", label = "quad")
ax2.scatter(lredshifts, la0, marker="o",color="darkviolet", alpha=0.8, label = "lin")
ax2.legend(loc=4)
ax2.set_ylabel("slope")
ax2.set_xlabel("redshift")
# plt.setp(ax2.get_xticklabels(), visible=False)
ax2.set_xlim([0.2,0.7])
plt.savefig("compare_lin_quad.pdf")


'''linslope, linintercept, r_value, p_value, std_err = linregress(lredshifts,lfactors0)
qslope, qintercept, r_value, p_value, std_err = linregress(qredshifts,qfactors0)

linfit=linslope*exes + linintercept
quadfit = qslope * exes + qintercept

plt.plot(exes, linfit, color = "purple")
plt.plot(exes, quadfit, color = "red")

annotation="""$linear$ $truth = {0:5.3f}*z +{1:5.3f}$
$quadratic$ $truth = {2:5.3f}*z +{3:5.3f}$
""".format(linslope, linintercept, qslope, qintercept)
plt.xlabel("Redshift")
plt.ylabel("Truth")
plt.title("Comparison of linear and quadratic models")

plt.text(0.2, 0.75, annotation, fontsize=12)
plt.legend(loc=4)

plt.show()
plt.savefig("compare_truths.pdf")'''
