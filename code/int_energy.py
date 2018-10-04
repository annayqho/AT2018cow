""" Compute the radiated energy (integrated over time)
and compare to store of energy computed on Day 22 """

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/annaho/Dropbox/Projects/Research/AT2018cow/code')
from astropy.cosmology import Planck15
from xray_lc import get_xrt, get_nustar
from scale_fluxes import sma_lc


d = Planck15.luminosity_distance(z=0.014).cgs.value

# Energy radiated in X rays
x,y,yerr = get_xrt() # in erg/cm^2/s

# Choose date range
choose = x < 22

# convert x to seconds
xsec = x[choose]*86400
xint = np.trapz(y[choose], xsec) # erg/cm2
xenrad = xint * 4 * np.pi * d**2
print(xenrad)

# 6.4E48 erg of energy before 20 days
# 3.4E48 erg radiated after 20 days

# Energy radiated at 230 GHz
a, b = sma_lc()
dt, f, ef = b
choose = dt < 22
dtsec = dt[choose]*86400
fcgs = f[choose] * 1E-3 * 1E-23
lcgs = 231.5E9 * fcgs
smaint = np.trapz(lcgs, dtsec) # erg/cm2 
smaenrad = smaint * 4 * np.pi * d**2
print(smaenrad)

# 5.1E46 erg radiated before 20 days
# 5.4E46 erg radiated after 20 days

