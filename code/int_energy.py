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

# convert x to seconds
xsec = x*86400
xint = np.trapz(y, xsec) # erg/cm2
xenrad = xint * 4 * np.pi * d**2

# 9.9E48 erg of energy as radiated in X rays

# Energy radiated at 230 GHz
a, b = sma_lc()
dt, f, ef = b
dtsec = dt*86400
fcgs = f * 1E-3 * 1E-23
lcgs = 231.5E9 * fcgs
smaint = np.trapz(lcgs, dtsec) # erg/cm2 
smaenrad = smaint * 4 * np.pi * d**2

# 1.2E47 erg of energy radiated at 231.5 GHz

# total SMA 231.5 and XRT: 1E49 erg
# and the minimum energy we get is 1.5E49 erg
