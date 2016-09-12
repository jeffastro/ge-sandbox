import numpy as np
from astropy.io import fits

data=fits.getdata('matchedkavli112315.fits')
columnnames=data.names

map_list = [(el, 1) for el in columnnames]

np.savetxt('matchedkavli112315_column_names.txt', map_list, fmt="%s")