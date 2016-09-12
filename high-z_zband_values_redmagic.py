from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt


filename = "stripe82_run_redmapper_v6.3.1_redmagic_1.0-02_wdr12.fit"
thecat = fits.getdata(filename)
redshift = thecat['ZREDMAGIC'].copy()
absmag = thecat['MODEL_MAG'].copy()
zband = absmag[:,3]

emceedata = fits.getdata("emceeRun031616_quad.fits")
intercepts = emceedata['a'].copy()
lins = emceedata['b'].copy()
quads = emceedata['c'].copy()
emceeredshifts = emceedata['redshift'].copy()

pivot = 19.8028

zbin = 0.64
halfwidth=0.01

zbandofinterest = zband[(redshift <= (zbin + halfwidth)) & (redshift >= (zbin - halfwidth))].flatten()

lowestzband = np.min(zbandofinterest)
inter = intercepts[emceeredshifts==0.64]
lin = lins[emceeredshifts==0.64]
quad = quads[emceeredshifts==0.64]

mstar = inter + lin*(lowestzband-pivot) + quad * (lowestzband-pivot) ** 2

print lowestzband
print mstar