#Modified 03.18.16 at 18:51

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt

matched_fits = fits.getdata("matched_cs82_redmag_mgc_031716.fits")
RA = matched_fits['DEC'].copy()
RAcs = matched_fits['DECcs'].copy()
plt.figure()
plt.scatter(RA, RAcs)
plt.xlabel("MGC/redMaGiC DEC")
plt.ylabel("CS82 DEC")
plt.show()