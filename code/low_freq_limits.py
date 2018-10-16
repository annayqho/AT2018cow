""" Get low-frequency VLASS limits """

import numpy as np
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.cosmology import Planck15
from get_radio import get_data_all
from synchrotron_fit import self_abs,fit_self_abs

# ultimately, we want a 3 GHz light curve
dt = []
lc = []

# Step 1: get all of the days with ATCA data
tel, freq, days, flux, eflux_form, eflux_sys = get_data_all()
ef = np.sqrt(eflux_form**2 + eflux_sys**2)
# there are 7 days with ATCA data
# of those 7 days, the first three have multiple points
# and the last four only have 34 GHz

# Step 2: for each of those days, fit nu squared and extrapolate
# down to 3 GHz, to generate a 3 GHz light curve
udays = np.unique(days[tel=='ATCA'])
for day in udays:
    dt.append(day)
    choose = np.logical_and(days == day, eflux_form != 99) # no upper limits
    popt, pcov = fit_self_abs(freq[choose], flux[choose], ef[choose])
    lc.append(self_abs(3, popt[0]))
dt = np.array(dt)
lc = np.array(lc)

# Step 3: Plot the 3 GHz light curve as luminosity
# give them all 20% uncertainties
dcm = Planck15.luminosity_distance(z=0.014).cgs.value
lum = lc*1E-3*1E-23*4*np.pi*dcm**2
plt.errorbar(dt, lum, yerr=0.2*lum, c='k', fmt='.')

# Step 4: Add in the VLASS limits

# Formatting
plt.tick_params(axis='both', labelsize=14)
plt.xscale('log')
plt.yscale('log')
plt.ylabel("$L_{\\nu}$ [erg/s/Hz]", fontsize=16)
plt.xlabel("Days since MJD 58285.0 (2018 June 16 UT)", fontsize=16)

plt.show()
