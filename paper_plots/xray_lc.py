""" Plot the X-ray light curve vs. the SMA light curves """

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from plot_spec import get_data
from plot_lc import sma

direc = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data/SwiftXRT/lc"
dat = Table.read(direc + "/table_to_plot.dat", format='ascii')
t = dat['!Time']
cts = dat['Rate']
ects =dat['Rate_Err']

# mean conversion factor
conv = 4.80e-11 # flux / cps
xrt_flux = cts * conv

fig, axarr = plt.subplots(2, 1, sharex=True)
ax = axarr[0]
days, nulow, nuhigh, flow, eflow, fhigh, efhigh = sma(ax, legend=False)

for ii in np.arange(14,19):
    ax.axvline(x=ii, alpha=0.5, lw=2)

# T0 for this burst is 2018 Jun 19 at 10:31:55.845 UT
# whereas I'm using June 16 as my T0
ax = axarr[1]
ax.errorbar(t/86400 + 3, cts, yerr=ects, fmt='.', c='k')
ax.set_xlabel("Days Since 16 June UT", fontsize=14)
ax.set_ylabel("XRT Counts/Sec", fontsize=14)

for ii in np.arange(14,19):
    ax.axvline(x=ii, alpha=0.5, lw=2)

#plt.show()
plt.savefig("xray_vs_radio.png")
