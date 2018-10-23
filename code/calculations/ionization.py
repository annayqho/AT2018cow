""" Calculate the ionization parameter

Let's day Day 22. We assume that it is a blackbody superimposed with the power
law. The blackbody is obviously dominant, so we can just add the two of them
together.
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/code")
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/code/paper_plots")
from get_optical import get_opt_sed,get_bb,get_nonthermal
from xray_lc import get_xrt, get_nustar
from astropy.modeling.blackbody import blackbody_nu
from astropy.cosmology import Planck15

# the range of frequencies over which we are doing our integral
# 13.6 eV = 3.3E15 keV
nu = np.logspace(np.log10(3.3E15), np.log10(2.4E19),10000)
d = Planck15.luminosity_distance(z=0.014).cgs.value

# constants
h = 4.14E-15
k = 1.38E-16
c = 3E10

# Blackbody
T = 16689
Lnu_bb = (2.9E30 * blackbody_nu(nu, T)).value
Fnu_bb = Lnu_bb / (4*np.pi*d**2)
#plt.plot(nu, nu*Lnu_bb, c='k')

# Power law
Lnu_pl = (1.6E24 * (nu/4.2E17)**(-0.54))
Fnu_pl = Lnu_pl / (4*np.pi*d**2)
#plt.plot(nu, nu*Lnu_pl, c='lightblue')

Lnu = Lnu_bb + Lnu_pl
Fnu = Fnu_bb + Fnu_pl
int_top = np.trapz(nu*Lnu, nu)
print(int_top)
int_bottom = np.trapz(Lnu, nu)
print(int_bottom)

temp = (int_top/int_bottom)/(4*k)

# plt.ylim(1E35, 1E44)
# plt.xlim(1E10, 1E19)
# plt.xscale('log')
# plt.yscale('log')
# plt.show()
