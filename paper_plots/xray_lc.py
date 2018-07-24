""" Plot the X-ray light curve vs. the SMA light curves """

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from plot_spec import get_data
from plot_lc import sma

direc = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
dat = Table.read(direc + "/xray_lc.txt", format='ascii', delimiter='&')
t = dat['t']
cts = dat['Lum']
ects =dat['elum']

fig, axarr = plt.subplots(2, 1, sharex=True)
ax = axarr[0]
days, nulow, nuhigh, flow, eflow, fhigh, efhigh = sma(ax, legend=False)

ax.text(
        0.1, 0.9, "SMA", 
        transform=ax.transAxes, fontsize=14, verticalalignment='top')

for ii in [12.9, 21.4, 24.4]:
    ax.axvline(x=ii, alpha=0.5, lw=2)

# T0 for this burst is 2018 Jun 19 at 10:31:55.845 UT
# whereas I'm using June 16 as my T0

# For the SMA data, the first date is
# 58290.00000000 + 0.391 # (9.393/24)
# 58290.391

ax = axarr[1]
ax.errorbar((t-t[0])+3, cts, yerr=ects, fmt='.', c='k')
ax.set_xlabel("Days Since 16 June UT", fontsize=14)
ax.set_ylabel(
r"$L (10^{43} \, \mathrm{erg s}^{-1} \, \mathrm{cm}^{-2})$", fontsize=14)

for ii in [12.9, 21.4, 24.4]:
    ax.axvline(x=ii, alpha=0.5, lw=2)

ax.text(
        0.1, 0.9, "Swift/XRT", 
        transform=ax.transAxes, fontsize=14, verticalalignment='top')

#plt.show()
plt.savefig("xray_vs_radio.png")
