# Last modified 03.27.16 at 19:35
# loop for creating many quadratic fits across z bins with a single pivot for all.
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# ------get data-------------------------------------------------------
thecat = fits.getdata("stripe82_run_redmapper_v6.3.1_redmagic_0.5-10_wdr12.fit")
thecat2 = fits.getdata("stripe82_run_redmapper_v6.3.1_redmagic_1.0-02_wdr12.fit")
emceedata = fits.getdata("emceeRun031616_quad.fits")
redshift = thecat['ZREDMAGIC'].copy()
model_mag = thecat['MODEL_MAG'].copy()
zband = model_mag[:,3]
redshift2 = thecat2['ZREDMAGIC'].copy()
model_mag2 = thecat2['MODEL_MAG'].copy()
zband2 = model_mag2[:,3]

intercepts = emceedata['a'].copy()
lins = emceedata['b'].copy()
quads = emceedata['c'].copy()
emceeredshifts = emceedata['redshift'].copy()

#get pivot point
kavli_fits = fits.getdata("matchedkavli112315.fits")
zbandforpivot = kavli_fits['ZBAND'].copy()
pivot = np.median(zbandforpivot)
print pivot

columnnames=['ID', 'RA', 'DEC', 'REDMAGICFLAG', 'IMAG', 'IMAG_ERR', 'MODEL_MAG', 'MODEL_MAGERR', 'MABS', 'MABS_ERR', 'ILUM', 'ZREDMAGIC', 'ZREDMAGIC_E', 'CHISQ', 'ZSPEC', 'PHOTOID']
columnformats=["J", "D", "D", "B", "E", "E", "4E", "4E", "4E", "4E", "E", "E", "E", "E", "E", "K"] 

mass1 = np.zeros(len(zband)).tolist()
mass2 = np.zeros(len(zband2)).tolist()

print "length of zband1: " + str(len(zband))
print "length of zband2: " + str(len(zband2))
for gal in xrange(len(zband)):
	currentz = redshift[gal]
	mcz = min(emceeredshifts, key=lambda x:abs(x-currentz))
	mass1[gal]=intercepts[emceeredshifts==mcz] + lins[emceeredshifts==mcz]*(zband[gal]-pivot) + quads[emceeredshifts==mcz]*(zband[gal]-pivot)**2


for gal in xrange(len(zband2)):
	currentz = redshift2[gal]
	mcz = min(emceeredshifts, key=lambda x:abs(x-currentz))
	mass2[gal]=intercepts[emceeredshifts==mcz] + lins[emceeredshifts==mcz]*(zband2[gal]-pivot) + quads[emceeredshifts==mcz]*(zband2[gal]-pivot)**2

cols=[]
for colnum in xrange(len(columnnames)):
	newarray=np.concatenate([thecat[columnnames[colnum]].copy(),thecat2[columnnames[colnum]].copy()], axis=0)
	cols.append(fits.Column(name=columnnames[colnum], format=columnformats[colnum], array=newarray))
	print newarray.shape

mass = mass1+mass2
mass = np.array(mass).flatten()
print "shape of mass: " + str(mass.shape)

cols.append(fits.Column(name='MCMC_MSTAR', format="E", array=mass))
cols = fits.ColDefs(cols)

tbhdu = fits.BinTableHDU.from_columns(cols)
tbhdu.writeto('stripe82_redmagic_with_stellarmass_sasha_both_cats.fits')
