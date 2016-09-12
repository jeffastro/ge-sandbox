import numpy as np
from astropy.io import fits

# ------get data-------------------------------------------------------
thecat = fits.getdata("matched_cs82_redmag_051916.fits")

headers = thecat.dtype.names
columns = []
redshift = thecat['ZREDMAGIC']

for header in headers:
    currentcolumn = thecat[header]
    galsofinterest = currentcolumn[(redshift<=(0.26)) & (redshift>=(0.24))] 
    columns.append( fits.Column(name=header, format = 'E', array = galsofinterest) )

tbhdu = fits.BinTableHDU.from_columns(columns)
tbhdu.writeto('middleredshiftgalaxies.fits')
