from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt

catpath = raw_input("catalog path and filename? ")
#mode = raw_input("linear or quadratic? ").lower()

#if (mode == "linear") or (mode == "lin") or (mode == "l"): linearmode = True
#else: linearmode = False

kavli_fits = fits.getdata(catpath)
mstar = kavli_fits['MASS_IR_BEST'].copy()
zband = kavli_fits['ZBAND'].copy()
redshift = kavli_fits['ZREDMAGIC'].copy()
# emass = kavli_fits['MASSERR_IR_BEST'].copy()
# ezband = kavli_fits['ZBANDERR'].copy()

#print np.max(redshift)
# while (zbin < limitontheoperation):    
bins = np.linspace(0.38, 0.6, num=10)
print
print bins
print
print
for zbin in bins:
    width = 0.02
    halfwidth = width/2.0
    medz_mstar = mstar[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth))]# & (mstar > 0) & (zband > 0)]
    medz_zband = zband[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth))]# & (mstar > 0) & (zband > 0)]
    # emasscurrent = emass[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (mstar > 0) & (zband > 0)]
    # ezbandcurrent = ezband[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth)) & (mstar > 0) & (zband > 0)]
    
    midzband = 20.0

    if len(medz_mstar)>1:
        #print medz_mstar
        bin1 = medz_mstar[(medz_zband<(midzband+0.1)) & (medz_zband>midzband) & (medz_mstar>0)]
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        plt.title("Redshift {0:1f} around z-band ".format(zbin) + str(midzband))
        ax1.set_xlabel("Stellar mass")
        ax1.set_ylabel("number of galaxies")
        ax1.hist(np.array(bin1), bins=10,align='left',color = "salmon")
        #plt.show()
        plt.savefig("choosing_estimator_z{0:2f}.png".format(zbin))
        plt.show()
    else:
        print "medz_mstar: " + str(medz_mstar)
        print "zbin: " + str(zbin)
