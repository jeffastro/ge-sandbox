#Last updated on 11.23.15 at 15:35
#This matches the KAVLI (Kevin Bundy) catalog to both Stripe82 catalogs

from astropy.coordinates import SkyCoord, Angle
from astropy.io import fits
from astropy import units as u
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


#filename1 = "/catalogs/kavli_catalog/UKBOSS_best_ukwide_v5_2-02jun2015.fits"
filename2 = "/calvin1/erozo/redmagic_catalogs/s82/v6.3.1/stripe82_run_redmapper_v6.3.1_redmagic_1.0-02_wdr12.fit"
filename3 = "/calvin1/erozo/redmagic_catalogs/s82/v6.3.1/stripe82_run_redmapper_v6.3.1_redmagic_0.5-10_wdr12.fit"
#filename1 = "/home/safonova/stcandle/matchingcats/kavli_lowerangles.fits"
filename1 = "/home/safonova/stcandle/kavli_catalog/UKBOSS_best_ukwide_v5_2-02jun2015.fits"

redmagic_fits = fits.getdata(filename2)
kavli_fits = fits.getdata(filename1)
redtwo_fits = fits.getdata(filename3)
print "Got the data"

RA = redmagic_fits['RA'].copy()
DEC = redmagic_fits['DEC'].copy()
ZREDMAGIC = redmagic_fits['ZREDMAGIC'].copy()
MODEL_MAG = redmagic_fits['MODEL_MAG'].copy()
ZBAND = MODEL_MAG[:,3].copy()
MODEL_MEGERR = redmagic_fits['MODEL_MAGERR'].copy()
ZBANDERR = MODEL_MAGERR[:,3]
ILUM = redmagic_fits['ILUM'].copy()

OBJID = kavli_fits['OBJID'].copy()
RAkavli = kavli_fits['RA'].copy()
DECkavli = kavli_fits['DEC'].copy()
ORIGIN_BEST = kavli_fits['ORIGIN_BEST'].copy()
ZBEST = kavli_fits['ZBEST'].copy()
ZPERT = kavli_fits['ZPERT'].copy()
MASS_IR_BEST = kavli_fits['MASS_IR_BEST'].copy()
MASSERR_IR_BEST = kavli_fits['MASSERR_IR_BEST'].copy()
MASS_IR_PERT = kavli_fits['MASS_IR_PERT'].copy()
MASS_OPT_BEST = kavli_fits['MASS_OPT_BEST'].copy()
MASSERR_OPT_BEST = kavli_fits['MASSERR_OPT_BEST'].copy()
MASS_OPT_PERT = kavli_fits['MASS_OPT_PERT'].copy()
B1000_IR_BEST = kavli_fits['B1000_IR_BEST'].copy()
B300_IR_BEST = kavli_fits['B300_IR_BEST'].copy()
KCOR_MASS_BEST = kavli_fits['KCOR_MASS_BEST'].copy()
ABSMAG_BEST = kavli_fits['ABSMAG_BEST'].copy()
SIGMAZ_REIS = kavli_fits['SIGMAZ_REIS'].copy()
SIGMAZ_RM = kavli_fits['SIGMAZ_RM'].copy()
SIGMAZ_ZRED = kavli_fits['SIGMAZ_ZRED'].copy()
SIGMAZ_BEST  = kavli_fits['SIGMAZ_BEST'].copy()
LUMERR_DEX = kavli_fits['LUMERR_DEX'].copy()


print "Turned the RA and DEC data into arrays for the first two catalogs"

red_ra = Angle(RA, unit=u.degree)                                              
red_dec = Angle(DEC,unit=u.degree)        
red_coord = SkyCoord(ra=red_ra,dec=red_dec, unit = (u.deg,u.deg))
print "shaped the first redmagic angles into red_coord"

kavli_ra = Angle(RAkavli, unit=u.degree)                               
kavli_dec = Angle(DECkavli,unit=u.degree)
kavli_coord = SkyCoord(ra=kavli_ra,dec=kavli_dec, unit = (u.deg,u.deg))
print "shaped the kavli angles into red_coord"

idx, d2d, d3d = red_coord.match_to_catalog_sky(kavli_coord)

print "idx is %s long" % (len(idx))

d2dvalues = d2d.value
print "d2d is %s long" %(len(d2d))
print "the mean d2d value is %s" %np.mean(d2d)
print "or %s" %np.mean(d2dvalues)
matchedRA = []
matchedDEC = []
matchedOBJID = []
matchedORIGIN = []
matchedZBEST = []
matchedZBEST = []
matchedZPERT = []
matchedMASS_IR_BEST = []
matchedMASSERR_IR_BEST = []
matchedMASS_IR_PERT = []
matchedMASS_OPT_BEST = []
matchedMASSERR_OPT_BEST = []
matchedMASS_OPT_PERT = []
matchedB1000_IR_BEST = []
matchedB300_IR_BEST = []
matchedKCOR_MASS_BEST = []
matchedABSMAG_BEST = []
matchedABSMAG_BEST = np.array(matchedABSMAG_BEST)
emptynine = [None, None, None, None, None, None, None, None, None]
emptynine = np.array(emptynine)
emptynine = emptynine.flatten()
matchedSIGMAZ_REIS = []
matchedSIGMAZ_RM = []
matchedSIGMAZ_ZRED = []
matchedSIGMAZ_BEST = []
matchedLUMERR_DEX = []        

for j in xrange(len(idx)):
    if d2dvalues[j] <= 0.000138888889:        
        matchedRA.append(RAkavli[idx[j]])
        matchedDEC.append(DECkavli[idx[j]])
#        matchedOBJID.append(OBJID[idx[j]])
        matchedORIGIN.append(ORIGIN_BEST[idx[j]])
        matchedZBEST.append(ZBEST[idx[j]])
        matchedZPERT.append(ZPERT[idx[j]])
        matchedMASS_IR_BEST.append(MASS_IR_BEST[idx[j]])
        matchedMASSERR_IR_BEST.append(MASSERR_IR_BEST[idx[j]])
        matchedMASS_IR_PERT.append(MASS_IR_PERT[idx[j]])
        matchedMASS_OPT_BEST.append(MASS_OPT_BEST[idx[j]])
        matchedMASSERR_OPT_BEST.append(MASSERR_OPT_BEST[idx[j]])
        matchedMASS_OPT_PERT.append(MASS_OPT_PERT[idx[j]])
        matchedB1000_IR_BEST.append(B1000_IR_BEST[idx[j]])
        matchedB300_IR_BEST.append(B300_IR_BEST[idx[j]])
        matchedKCOR_MASS_BEST.append(KCOR_MASS_BEST[idx[j]])
        np.concatenate((matchedABSMAG_BEST, ABSMAG_BEST[idx[j]]), axis=0)
        matchedSIGMAZ_REIS.append(SIGMAZ_REIS[idx[j]])
        matchedSIGMAZ_RM.append(SIGMAZ_RM[idx[j]])
        matchedSIGMAZ_ZRED.append(SIGMAZ_ZRED[idx[j]])
        matchedSIGMAZ_BEST.append(SIGMAZ_BEST[idx[j]])
        matchedLUMERR_DEX.append(LUMERR_DEX[idx[j]])
    else:
        matchedRA.append(None)
        matchedDEC.append(None)
#        matchedOBJID.append(None)
        matchedORIGIN.append(None)
        matchedZBEST.append(None)
        matchedZPERT.append(None)
        matchedMASS_IR_BEST.append(None)
        matchedMASSERR_IR_BEST.append(None)
        matchedMASS_IR_PERT.append(None)
        matchedMASS_OPT_BEST.append(None)
        matchedMASSERR_OPT_BEST.append(None)
        matchedMASS_OPT_PERT.append(None)
        matchedB1000_IR_BEST.append(None)
        matchedB300_IR_BEST.append(None)
        matchedKCOR_MASS_BEST.append(None)
        np.concatenate((matchedABSMAG_BEST, emptynine),axis=0)
        matchedSIGMAZ_REIS.append(None)
        matchedSIGMAZ_RM.append(None)
        matchedSIGMAZ_ZRED.append(None)
        matchedSIGMAZ_BEST.append(None)
        matchedLUMERR_DEX.append(None)
print("matched the kavli to the first redmagic catalog")
print "matchedORIGIN is now %s long" %len(matchedORIGIN)
print "the first element in idx is:"
print idx[0]

#-------------------------Second catalog-------------------------

RA2 = redtwo_fits['RA'].copy()
RA2 = np.array(RA2)

DEC2 = redtwo_fits['DEC'].copy()
DEC2 = np.array(DEC2)

ZREDMAGIC2 = redtwo_fits['ZREDMAGIC'].copy()
ZREDMAGIC2 = np.array(ZREDMAGIC2)
ILUM2 = redtwo_fits['ILUM'].copy()
ILUM2 = np.array(ILUM2)
MODEL_MAG2 = redtwo_fits['MODEL_MAG'].copy()
ZBAND2 = MODEL_MAG2[:,3]
MODEL_MEGERR2 = redtwo_fits['MODEL_MAGERR'].copy()
ZBANDERR2 = MODEL_MAGERR2[:,3]

redtwo_ra = Angle(RA2, unit=u.degree)                                              
redtwo_dec = Angle(DEC2,unit=u.degree)        
redtwo_coord = SkyCoord(ra=redtwo_ra,dec=redtwo_dec, unit = (u.deg,u.deg))
print "shaped the second redmagic angles into red_coord"

idx, d2d, d3d = redtwo_coord.match_to_catalog_sky(kavli_coord)

print "this idx is %s long" % (len(idx))

d2dvalues = d2d.value
print "the mean d2d value is %s" %np.mean(d2d)
print "or %s" %np.mean(d2dvalues)
print "the min d2d value is %s" %np.min(d2dvalues)
print "the median d2d value is %s" %np.median(d2dvalues)

for j in xrange(len(idx)):
    if d2dvalues[j] <= 0.000138888889:
#        matchedOBJID.append(OBJID[idx[j]])
        matchedRA.append(RAkavli[idx[j]])
        matchedDEC.append(DECkavli[idx[j]])
        matchedORIGIN.append(ORIGIN_BEST[idx[j]])
        matchedZBEST.append(ZBEST[idx[j]])
        matchedZPERT.append(ZPERT[idx[j]])
        matchedMASS_IR_BEST.append(MASS_IR_BEST[idx[j]])
        matchedMASSERR_IR_BEST.append(MASSERR_IR_BEST[idx[j]])
        matchedMASS_IR_PERT.append(MASS_IR_PERT[idx[j]])
        matchedMASS_OPT_BEST.append(MASS_OPT_BEST[idx[j]])
        matchedMASSERR_OPT_BEST.append(MASSERR_OPT_BEST[idx[j]])
        matchedMASS_OPT_PERT.append(MASS_OPT_PERT[idx[j]])
        matchedB1000_IR_BEST.append(B1000_IR_BEST[idx[j]])
        matchedB300_IR_BEST.append(B300_IR_BEST[idx[j]])
        matchedKCOR_MASS_BEST.append(KCOR_MASS_BEST[idx[j]])
        np.concatenate((matchedABSMAG_BEST, ABSMAG_BEST[idx[j]]), axis=0)
        matchedSIGMAZ_REIS.append(SIGMAZ_REIS[idx[j]])
        matchedSIGMAZ_RM.append(SIGMAZ_RM[idx[j]])
        matchedSIGMAZ_ZRED.append(SIGMAZ_ZRED[idx[j]])
        matchedSIGMAZ_BEST.append(SIGMAZ_BEST[idx[j]])
        matchedLUMERR_DEX.append(LUMERR_DEX[idx[j]])
    else:
#        matchedOBJID.append(None)
        matchedRA.append(None)
        matchedDEC.append(None)
        matchedORIGIN.append(None)
        matchedZBEST.append(None)
        matchedZPERT.append(None)
        matchedMASS_IR_BEST.append(None)
        matchedMASSERR_IR_BEST.append(None)
        matchedMASS_IR_PERT.append(None)
        matchedMASS_OPT_BEST.append(None)
        matchedMASSERR_OPT_BEST.append(None)
        matchedMASS_OPT_PERT.append(None)
        matchedB1000_IR_BEST.append(None)
        matchedB300_IR_BEST.append(None)
        matchedKCOR_MASS_BEST.append(None)
        np.concatenate((matchedABSMAG_BEST, emptynine), axis=0)
        matchedSIGMAZ_REIS.append(None)
        matchedSIGMAZ_RM.append(None)
        matchedSIGMAZ_ZRED.append(None)
        matchedSIGMAZ_BEST.append(None)
        matchedLUMERR_DEX.append(None)
        
print("matched the second redmagic catalog to kavli")
print("matchedORIGIN is %s long") %len(matchedORIGIN)
print "d2d is %s long" %len(d2d)
print "the first element in d2d is %s" % d2d[0]

RA = np.insert(RA2, np.arange(len(RA)), RA)
DEC = np.insert(DEC2, np.arange(len(DEC)), DEC)
ZREDMAGIC = np.insert(ZREDMAGIC2, np.arange(len(ZREDMAGIC)), ZREDMAGIC)
ILUM = np.insert(ILUM2, np.arange(len(ILUM)), ILUM)
ZBAND = np.insert(ZBAND2, np.arange(len(ZBAND)), ZBAND)
ZBANDERR = np.insert(ZBANDERR2, np.arange(len(ZBANDERR)), ZBANDERR)

print("connected the two redmagic catalogs together")

#c1=fits.Column(name = 'OBJID', format = 'I', array = np.array(matchedOBJID))
c2=fits.Column(name = 'RAkavli', format = 'E', array = np.array(matchedRA))
c3=fits.Column(name = 'DECkavli', format = 'E', array = np.array(matchedDEC))
c4=fits.Column(name = 'ORIGIN_BEST', format = 'E', array = np.array(matchedORIGIN))
c5=fits.Column(name = 'ZBESTkavli', format = 'E', array = np.array(matchedZBEST))
c6=fits.Column(name = 'ZPERT', format = 'E', array = np.array(matchedZPERT))
c7=fits.Column(name = 'MASS_IR_BEST', format = 'E', array = np.array(matchedMASS_IR_BEST))
c8=fits.Column(name = 'MASSERR_IR_BEST', format = 'E', array = np.array(matchedMASSERR_IR_BEST))
c9=fits.Column(name = 'MASS_IR_PERT', format = 'E', array = np.array(matchedMASS_IR_PERT))
c10=fits.Column(name = 'MASS_OPT_BEST', format = 'E', array = np.array(matchedMASS_OPT_BEST))
c11=fits.Column(name = 'MASSERR_OPT_BEST', format = 'E', array = np.array(matchedMASSERR_OPT_BEST))
c12=fits.Column(name = 'MASS_OPT_PERT', format = 'E', array = np.array(matchedMASS_OPT_PERT))
c13=fits.Column(name = 'B1000_IR_BEST', format = 'E', array = np.array(matchedB1000_IR_BEST))
c14=fits.Column(name = 'B300_IR_BEST', format = 'E', array = np.array(matchedB300_IR_BEST))
c15=fits.Column(name = 'KCOR_MASS_BEST', format = 'E', array = np.array(matchedKCOR_MASS_BEST))
c16=fits.Column(name = 'ABSMAG_BEST', format = '9E', array = np.array(matchedABSMAG_BEST))
c17=fits.Column(name = 'SIGMAZ_REIS', format = 'E', array = np.array(matchedSIGMAZ_REIS))
c18=fits.Column(name = 'SIGMAZ_RM', format = 'E', array = np.array(matchedSIGMAZ_RM))
c19=fits.Column(name = 'SIGMAZ_ZRED', format = 'E', array = np.array(matchedSIGMAZ_ZRED))
c20=fits.Column(name = 'SIGMAZ_BEST', format = 'E', array = np.array(matchedSIGMAZ_BEST))
c21=fits.Column(name = 'LUMERR_DEX', format = 'E', array = np.array(matchedLUMERR_DEX))
c22=fits.Column(name = 'RA', format='E', array = np.array(RA))
c23=fits.Column(name = 'DEC', format='E', array = np.array(DEC))
c24=fits.Column(name = 'ZREDMAGIC', format='E', array = np.array(ZREDMAGIC))
c25=fits.Column(name = 'ZBAND', format='E', array = np.array(ZBAND))
c26=fits.Column(name = 'ZBANDERR', format='E', array = np.array(ZBANDERR))

tbhdu = fits.BinTableHDU.from_columns([c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26])
tbhdu.writeto('matchedkavli051716.fits')
