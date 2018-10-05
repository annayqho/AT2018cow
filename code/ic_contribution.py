""" Check for IC compton contribution at early times """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from astropy.cosmology import Planck15
from get_optical import *
from xray_lc import *
from scale_fluxes import sma_lc


dt,lum = get_opt_lum()

f = 0.5 # filling factor
t0 = 22
R0 = 7E15
c = 3E10
t = np.arange(1,30)
R = R0 * (t/t0)
V = (4/3) * f * np.pi * R**3
v = 0.13*c


# Radiation energy density
L = np.interp(t,dt,lum)
Uph = L / (4 * np.pi * R**2 * c)

# Magnetic energy density
UB = (6**2 / (8*np.pi))

# Ratio of energy densities
Uratio = Uph/UB

# Ratio of X-ray to radio luminosity
d = Planck15.luminosity_distance(z=0.014).cgs.value
dt_x, f_x, ef_x = get_xrt()
l_x = np.interp(t, dt_x, f_x * 4 * np.pi * d**2)

a, b = sma_lc()
dt_rad, f_rad, ef = b
ef_rad = np.sqrt(ef**2 + (0.15*f)**2)
l_rad = np.interp(t, dt_rad, 230E9 * f_rad * 1E-3 * 1E-23 * 4 * np.pi * d**2 )
Lratio = l_x/l_rad
