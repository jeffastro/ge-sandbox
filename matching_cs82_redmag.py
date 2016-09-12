#/calvin1/erozo/redmagic_catalogs/s82/v6.3.1/stripe82_run_redmapper_v6.3.1_redmagic_1.0-02_wdr12.fit#Matching the CS82 to the matched redmagic/MGC catalog
#Modified 04.03.16 at 15:39

from astropy.coordinates import SkyCoord, Angle
from astropy.io import fits
from astropy import units as u
import numpy as np

filename1="/calvin1/erozo/redmagic_catalogs/s82/v6.3.1/stripe82_run_redmapper_v6.3.1_redmagic_1.0-02_wdr12.fit"
filename2 = "/calvin1/erozo/redmagic_catalogs/CS82/cs82_morphology_ser.fit"
filename3 = "/calvin1/erozo/redmagic_catalogs/s82/v6.3.1/stripe82_run_redmapper_v6.3.1_redmagic_0.5-10_wdr12.fit"

redmag1_fits = fits.getdata(filename1)
cs82_fits = fits.getdata(filename2)
redmag2_fits = fits.getdata(filename3)
print "Got the data"

def getreddata(red_catalog_fits):
    RA = red_catalog_fits['RA'].copy()
    DEC = red_catalog_fits['DEC'].copy()
    ZREDMAGIC = red_catalog_fits['ZREDMAGIC'].copy()
    MODEL_MAG = red_catalog_fits['MODEL_MAG'].copy()
    ZBAND = MODEL_MAG[:,3]
    return RA, DEC, ZREDMAGIC, ZBAND

RAcs = cs82_fits['ALPHA_J2000'].copy()
DECcs = cs82_fits['DELTA_J2000'].copy()
MAG_AUTO = cs82_fits['MAG_AUTO'].copy()
MAGERR_AUTO = cs82_fits['MAGERR_AUTO'].copy()
FLUX_RADIUS = cs82_fits['FLUX_RADIUS'].copy()
FLUX_MODEL = cs82_fits['FLUX_MODEL'].copy()
SPHEROID_REFF_WORLD = cs82_fits['SPHEROID_REFF_WORLD'].copy()
SPHEROID_REFFERR_WORLD = cs82_fits['SPHEROID_REFFERR_WORLD'].copy()
SPHEROID_SERSICN = cs82_fits['SPHEROID_SERSICN'].copy()

def matchingroutine (RAred, DECred, RAcs, DECcs):
    redmag_ra = Angle(RAred, unit=u.degree)                                              
    redmag_dec = Angle(DECred,unit=u.degree)        
    redmag_coord = SkyCoord(ra=redmag_ra,dec=redmag_dec, unit = (u.deg,u.deg))
    print "shaped the matched angles into matched_coord"

    cs82_ra = Angle(RAcs, unit=u.degree)                               
    cs82_dec = Angle(DECcs,unit=u.degree)
    cs82_coord = SkyCoord(ra=cs82_ra,dec=cs82_dec, unit = (u.deg,u.deg))
    print "shaped the cs82 angles into cs82_coord"

    idx, d2d, d3d = redmag_coord.match_to_catalog_sky(cs82_coord)

    print "idx is %s long" % (len(idx))
    return idx, d2d.value


RA, DEC, ZREDMAGIC, ZBAND = getreddata(redmag1_fits)
print "Turned the RA and DEC data into arrays for the first two catalogs"

idx, d2dvalues = matchingroutine(RA, DEC, RAcs, DECcs)

print "or %s" %np.mean(d2dvalues)

matchedMAG_AUTO = []
matchedMAGERR_AUTO = []
matchedFLUX_RADIUS = []
matchedFLUX_MODEL = []
matchedSPHEROID_REFF_WORLD = []
matchedSPHEROID_REFFERR_WORLD = []
matchedSPHEROID_SERSICN = []
newRA = []
newDEC = []
newZREDMAGIC = []
newZBAND = []

for j in xrange(len(idx)):
    if d2dvalues[j] <= 0.000138888889:        
        matchedMAG_AUTO.append(MAG_AUTO[idx[j]])
        matchedMAGERR_AUTO.append(MAGERR_AUTO[idx[j]])
        matchedFLUX_RADIUS.append(FLUX_RADIUS[idx[j]])
        matchedFLUX_MODEL.append(FLUX_MODEL[idx[j]])
        matchedSPHEROID_REFF_WORLD.append(SPHEROID_REFF_WORLD[idx[j]])
        matchedSPHEROID_REFFERR_WORLD.append(SPHEROID_REFFERR_WORLD[idx[j]])
        matchedSPHEROID_SERSICN.append(SPHEROID_SERSICN[idx[j]])
        newRA.append(RA[j])
        newDEC.append(DEC[j])
        newZREDMAGIC.append(ZREDMAGIC[j])
        newZBAND.append(ZBAND[j])


RA2, DEC2, ZREDMAGIC2, ZBAND2 = getreddata(redmag2_fits)
idx, d2dvalues = matchingroutine(RA2, DEC2, RAcs, DECcs)
print "Turned the RA and DEC data into arrays and matched second two catalogs"

for j in xrange(len(idx)):
    if d2dvalues[j] <= 0.000138888889:        
        matchedMAG_AUTO.append(MAG_AUTO[idx[j]])
        matchedMAGERR_AUTO.append(MAGERR_AUTO[idx[j]])
        matchedFLUX_RADIUS.append(FLUX_RADIUS[idx[j]])
        matchedFLUX_MODEL.append(FLUX_MODEL[idx[j]])
        matchedSPHEROID_REFF_WORLD.append(SPHEROID_REFF_WORLD[idx[j]])
        matchedSPHEROID_REFFERR_WORLD.append(SPHEROID_REFFERR_WORLD[idx[j]])
        matchedSPHEROID_SERSICN.append(SPHEROID_SERSICN[idx[j]])
        newRA.append(RA2[j])
        newDEC.append(DEC2[j])
        newZREDMAGIC.append(ZREDMAGIC2[j])
        newZBAND.append(ZBAND2[j])

c1=fits.Column(name = 'ZREDMAGIC', format = 'E', array = np.array(newZREDMAGIC))
c2=fits.Column(name = 'ZBAND', format = 'E', array = np.array(newZBAND))
c3=fits.Column(name = 'MAG_AUTO', format = 'E', array = np.array(matchedMAG_AUTO))
c4=fits.Column(name = 'MAGERR_AUTO', format = 'E', array = np.array(matchedMAGERR_AUTO))
c5=fits.Column(name = 'FLUX_RADIUS', format = 'E', array = np.array(matchedFLUX_RADIUS))
c6=fits.Column(name = 'FLUX_MODEL', format = 'E', array = np.array(matchedFLUX_MODEL))
c7=fits.Column(name = 'SPHEROID_REFF_WORLD', format = 'E', array = np.array(matchedSPHEROID_REFF_WORLD))
c8=fits.Column(name = 'SPHEROID_REFFERR_WORLD', format = 'E', array = np.array(matchedSPHEROID_REFFERR_WORLD))
c9=fits.Column(name = 'SPHEROID_SERSICN', format = 'E', array = np.array(matchedSPHEROID_SERSICN))
c10=fits.Column(name = 'RA', format = 'E', array = np.array(newRA))
c11=fits.Column(name = 'DEC', format = 'E', array = np.array(newDEC))

 
tbhdu = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11])
tbhdu.writeto('matched_cs82_redmag_040316.fits')

