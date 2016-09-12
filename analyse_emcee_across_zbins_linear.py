from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.ticker

filename="emceeRun022816_linear.fits"

emceedata=fits.getdata(filename)
columnnames=emceedata.names

a0= emceedata['a'].copy()
b0 = emceedata['b'].copy()
sigmas0 = emceedata['sigma'].copy()
factors0 = emceedata['factor'].copy()
redshifts = emceedata['redshift'].copy()
galnumbers = emceedata['number_of_galaxies'].copy()

sigmas1 = emceedata['sigmalow'].copy()
sigmas2 = emceedata['sigmahigh'].copy()
a1= emceedata['a_low'].copy()
b1 = emceedata['b_low'].copy()
a2= emceedata['a_high'].copy()
b2 = emceedata['b_high'].copy()
factors1 = emceedata['factor_low'].copy()
factors2 = emceedata['factor_high'].copy()

'''plt.figure(figsize=(14, 9), dpi=80)
plt.title("Linear fit")
plt.tight_layout()
plt.style.use('ggplot')
ax1 = plt.subplot2grid((3,2), (0,0))
ax2 = plt.subplot2grid((3,2), (0,1))
ax3 = plt.subplot2grid((3,2), (1, 0))
ax4 = plt.subplot2grid((3,2), (1, 1))
ax5 = plt.subplot2grid((3,2), (2, 0))
ax6 = plt.subplot2grid((3,2), (2, 1))

plt.subplots_adjust(left=0.1, bottom=0.08, right=None, top=None,
                wspace=None, hspace=None)
plt.suptitle('Linear fit', fontsize=12)

ax2.errorbar(redshifts, a0, yerr=[a1,a2], linestyle="None", marker=".",markerfacecolor="skyblue", color="green")
#ax2.title("Slope vs redshift in a quadratic model of MSTAR vs ZBAND")
#ax2.set_xlabel("Redshift")
ax2.set_ylabel("a")
plt.setp(ax2.get_xticklabels(), visible=False)
#plt.setp(ax2.get_yticklabels(), fontsize=6 )
#ax1.savefig("slopevsredshift.pdf")


ax3.errorbar(redshifts, b0, yerr=[b1,b2], linestyle="None", marker=".",markerfacecolor="skyblue", color="green")
#ax3.title("b vs redshift")
#ax3.set_xlabel("Redshift")
ax3.set_ylabel("b")
plt.setp(ax3.get_xticklabels(), visible=False)
#plt.setp(ax3.get_yticklabels(), fontsize=6 )

ax4.errorbar(redshifts, c0, yerr=[c1,c2], linestyle="None", marker=".",markerfacecolor="skyblue", color="green")
#ax4.title("c vs redshift")
#ax4.set_xlabel("Redshift")
plt.setp(ax4.get_xticklabels(), visible=False)
#plt.setp(ax4.get_yticklabels(), fontsize=6 )
ax4.set_ylabel("c")

#ax5.figure()
ax5.errorbar(redshifts, sigmas0, yerr=[sigmas1,sigmas2], linestyle="None", marker=".",markerfacecolor="skyblue", color="green")
#ax5.title("Sigma vs redshift in a linear model of MSTAR vs ZBAND")
ax5.set_xlabel("Redshift")
ax5.set_ylabel("$\sigma$")
#plt.setp(ax5.get_yticklabels(), fontsize=6 )
#plt.setp(ax5.get_xticklabels(), fontsize=6 )


ax1.errorbar(redshifts, factors0, yerr=[factors1,factors2], linestyle="None", marker=".",markerfacecolor="skyblue", color="green")
ax1.set_ylabel("Truth factor")
plt.setp(ax1.get_xticklabels(), visible=False)
#plt.setp(ax1.get_yticklabels(), fontsize=10 )


ax6.scatter(redshifts, galnumbers, marker=".",color="green")
ax6.set_xlabel("Redshift")
ax6.set_ylabel("Galaxies")
ax6.set_xlim([0.2,0.7])
plt.setp(ax6.get_yticklabels(), fontsize=6 )
#plt.setp(ax6.get_xticklabels(), fontsize=10 )

#y_formatter = matplotlib.ticker.set_scientific(True)
#ax6.yaxis.set_major_formatter(y_formatter)
plt.savefig("emcee_quad_allbins.pdf")


plt.show()'''


interceptbestfit= np.polyfit(redshifts, b0, 1)
lowlin,lowinter=np.poly1d(interceptbestfit)
interceptbestfit = lowinter + lowlin*redshifts

'''plt.figure()
plt.errorbar(redshifts, b0, yerr=[a1,a2], linestyle="None", marker=".",markerfacecolor="skyblue", color="green")
plt.plot(redshifts, interceptbestfit)
plt.ylabel("intercept")
plt.xlabel("redshift")
# plt.setp(ax2.get_xticklabels(), visible=False)
plt.title("Linear mcmc intercepts\nslope of fit line = {0}".format(lowlin))
plt.show()'''

plt.figure()
plt.errorbar(redshifts, a0, yerr=[b1,b2], linestyle="None", marker=".",markerfacecolor="skyblue", color="green")
plt.plot(redshifts, np.mean(a0))
plt.ylabel("slope")
plt.xlabel("redshift")
plt.title("Linear mcmc slopes\nmean slope = {0}".format(np.mean(a0)))
plt.xlabel
plt.show()