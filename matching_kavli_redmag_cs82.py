#Matching the CS82 to the matched redmagic/MGC catalog
#Modified 03.18.16 at 18:05

from astropy.coordinates import SkyCoord, Angle
from astropy.io import fits
from astropy import units as u
import numpy as np
#import matplotlib as mpl
import matplotlib.pyplot as plt


filename1 = "/home/safonova/stcandle/matchingcats/matchedkavli112315.fits"
filename2 = "/calvin1/erozo/redmagic_catalogs/CS82/cs82_morphology_ser.fit"

matched_fits = fits.getdata(filename1)
cs82_fits = fits.getdata(filename2)
print "Got the data"


RA = matched_fits['RA'].copy()
DEC = matched_fits['DEC'].copy()
ZREDMAGIC = matched_fits['ZREDMAGIC'].copy()
MASS_IR_BEST = matched_fits['MASS_IR_BEST'].copy()
MASSERR_IR_BEST = matched_fits['MASSERR_IR_BEST'].copy()
ABSMAG_BEST = matched_fits['ABSMAG_BEST'].copy()
ZBAND = matched_fits['ZBAND'].copy()

RAcs = cs82_fits['ALPHA_J2000'].copy()
DECcs = cs82_fits['DELTA_J2000'].copy()
MAG_AUTO = cs82_fits['MAG_AUTO'].copy()
MAGERR_AUTO = cs82_fits['MAGERR_AUTO'].copy()
FLUX_RADIUS = cs82_fits['FLUX_RADIUS'].copy()
FLUX_MODEL = cs82_fits['FLUX_MODEL'].copy()
SPHEROID_REFF_WORLD = cs82_fits['SPHEROID_REFF_WORLD'].copy()
SPHEROID_REFFERR_WORLD = cs82_fits['SPHEROID_REFFERR_WORLD'].copy()
SPHEROID_SERSICN = cs82_fits['SPHEROID_SERSICN'].copy()

print "Turned the RA and DEC data into arrays for the first two catalogs"

matched_ra = Angle(RA, unit=u.degree)                                              
matched_dec = Angle(DEC,unit=u.degree)        
matched_coord = SkyCoord(ra=matched_ra,dec=matched_dec, unit = (u.deg,u.deg))
print "shaped the matched angles into matched_coord"

cs82_ra = Angle(RAcs, unit=u.degree)                               
cs82_dec = Angle(DECcs,unit=u.degree)
cs82_coord = SkyCoord(ra=cs82_ra,dec=cs82_dec, unit = (u.deg,u.deg))
print "shaped the cs82 angles into cs82_coord"

idx, d2d, d3d = matched_coord.match_to_catalog_sky(cs82_coord)

print "idx is %s long" % (len(idx))

d2dvalues = d2d.value
print "d2d is %s long" %(len(d2d))
print "the mean d2d value is %s" %np.mean(d2d)
print "or %s" %np.mean(d2dvalues)

matchedMAG_AUTO = []
matchedMAGERR_AUTO = []
matchedFLUX_RADIUS = []
matchedFLUX_MODEL = []
matchedSPHEROID_REFF_WORLD = []
matchedSPHEROID_REFFERR_WORLD = []
matchedSPHEROID_SERSICN = []
newRAcs = []
newDECcs = []

for j in xrange(len(idx)):
    if d2dvalues[j] <= 0.000138888889:        
        matchedMAG_AUTO.append(MAG_AUTO[idx[j]])
        matchedMAGERR_AUTO.append(MAGERR_AUTO[idx[j]])
        matchedFLUX_RADIUS.append(FLUX_RADIUS[idx[j]])
        matchedFLUX_MODEL.append(FLUX_MODEL[idx[j]])
        matchedSPHEROID_REFF_WORLD.append(SPHEROID_REFF_WORLD[idx[j]])
        matchedSPHEROID_REFFERR_WORLD.append(SPHEROID_REFFERR_WORLD[idx[j]])
        matchedSPHEROID_SERSICN.append(SPHEROID_SERSICN[idx[j]])
        newRAcs.append(RAcs[idx[j]])
        newDECcs.append(DECcs[idx[j]])
    else:
        matchedMAG_AUTO.append(None)
        matchedMAGERR_AUTO.append(None)
        matchedFLUX_RADIUS.append(None)
        matchedFLUX_MODEL.append(None)
        matchedSPHEROID_REFF_WORLD.append(None)
        matchedSPHEROID_REFFERR_WORLD.append(None)
        matchedSPHEROID_SERSICN.append(None)
        newRAcs.append(None)
        newDECcs.append(None)


print("matched the kavli to the first redmagic catalog")
print "the first element in idx is:"
print idx[0]

c1=fits.Column(name = 'RA', format = 'E', array = np.array(RA))
c2=fits.Column(name = 'DEC', format = 'E', array = np.array(DEC))
c3=fits.Column(name = 'ZREDMAGIC', format = 'E', array = np.array(ZREDMAGIC))
c4=fits.Column(name = 'MASS_IR_BEST', format = 'E', array = np.array(MASS_IR_BEST))
c5=fits.Column(name = 'MASSERR_IR_BEST', format = 'E', array = np.array(MASSERR_IR_BEST))
c6=fits.Column(name = 'ABSMAG_BEST ', format = '9E', array = np.array(ABSMAG_BEST ))
c7=fits.Column(name = 'ZBAND', format = 'E', array = np.array(ZBAND))
c8=fits.Column(name = 'MAG_AUTO', format = 'E', array = np.array(matchedMAG_AUTO))
c9=fits.Column(name = 'MAGERR_AUTO', format = 'E', array = np.array(matchedMAGERR_AUTO))
c10=fits.Column(name = 'FLUX_RADIUS', format = 'E', array = np.array(matchedFLUX_RADIUS))
c11=fits.Column(name = 'FLUX_MODEL', format = 'E', array = np.array(matchedFLUX_MODEL))
c12=fits.Column(name = 'SPHEROID_REFF_WORLD', format = 'E', array = np.array(matchedSPHEROID_REFF_WORLD))
c13=fits.Column(name = 'SPHEROID_REFFERR_WORLD', format = 'E', array = np.array(matchedSPHEROID_REFFERR_WORLD))
c14=fits.Column(name = 'SPHEROID_SERSICN', format = 'E', array = np.array(matchedSPHEROID_SERSICN))
c15=fits.Column(name = 'RAcs', format = 'E', array = np.array(newRAcs))
c16=fits.Column(name = 'DECcs', format = 'E', array = np.array(newDECcs))

 
tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16])
tbhdu.writeto('matched_cs82_redmag_mgc_031716.fits')
