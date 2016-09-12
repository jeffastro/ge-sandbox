# -*- coding: utf-8 -*-
"""Author: Sasha Safonova

FSPS code is based on:
* Conroy, Gunn, & White 2009, ApJ, 699, 486
* Conroy & Gunn 2010, ApJ, 712, 833

FSPS bindings provided by:
* Dan Foreman-Mackey, Jonathan Sick, Queen's University, Ben Johnson \
https://github.com/dfm/python-fsps
"""

import numpy as np
from matplotlib import pyplot as plt
import fsps

# sp = fsps.StellarPopulation(compute_vega_mags=False, fburst=1, dust2=0,
#                            tburst=5, sfh=4, add_stellar_remnants=1, tau=5,
#                            zred=0.4)
                            
sp = fsps.StellarPopulation(compute_vega_mags=False, fburst=1, dust2=0,
                            tburst=2.1, sfh=0, add_stellar_remnants=1)

masses = sp.stellar_mass
sdss_bands = fsps.find_filter('sdss')

mags = sp.get_mags(tage=13.7, bands=sdss_bands)
wave, spec = sp.get_spectrum(tage=13.7)

use_vega_mags, apply_redshift = 0, 0

imf, imf1, imf2, imf3, vdmc, mdave, dell, delt, sbss, fbhb, pagb = \
    0,1.3,2.3,2.3,0.08,0.5,0.,0.,0.,0.,1.
        
model_name = "test_fsps"
dust_type = 0
zmet = 0
sfh = 1
tau = 0.5
const = 0.1
fburst = 0.
tburst = 0.
dust_tesc = 7.
dust1 = 0.
dust2 = 0.
dust_clumps = -99.
frac_no_dust = 0.1
dust_index = -0.7,
mwr = 3.1
wgp1 = 1
wgp2 = 1
wgp3 = 1
duste_gamma = 0.01
duste_umin = 1.0
duste_qpah = 3.5
tage = 0.
        
'''fsps.driver.setup(use_vega_mags, apply_redshift)
nBands = fsps.driver.get_n_bands()
nAges = fsps.driver.get_n_ages()
nMasses = fsps.driver.get_n_masses()
# print "There are", nBands, "bands"
# print "There are", nAges, "ages"
print ("There are" + nMasses + "masses")'''