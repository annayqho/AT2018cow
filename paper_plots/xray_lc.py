""" Plot the X-ray light curve vs. the SMA light curves

X ray: 2018 Jun 19 at 10:34:30.742 UT was T0
First SMA track was 2018 Jun 21, 03:58:25 - 14:18:09 UTC (10h 19m)
I call this Day 5. That means that the XRT point is Day 3.
So I need to add 3 days to all XRT points.
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from plot_spec import get_data
from plot_lc import sma
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)

direc = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
dat = Table.read(direc + "/xray_lc.txt", format='ascii')
t = dat['col1'] / (3600*24) # in days
# count to flux conversion (absorbed): 4.26 × 10-11 erg cm-2 ct-1
cts = dat['col4'] * 4.26
ects =dat['col5']

fig, axarr = plt.subplots(2, 1, sharex=True)
ax = axarr[0]
dat = Table.read(
    "../data/radio_lc.dat", delimiter="&", format='ascii.no_header')
tel = np.array(dat['col2'])
choose = tel == 'SMA'

days = np.array(dat['col1'][choose])
freq = np.array(dat['col3'][choose]).astype(float)
flux_raw = np.array(dat['col4'][choose])
flux = np.array(
        [float(val.split("pm")[0][1:]) for val in flux_raw])
eflux = np.array(
        [float(val.split("pm")[1][0:-1]) for val in flux_raw])

choose = np.logical_and(freq > 230, freq < 245)
low_freq = min(freq[choose])
max_freq = max(freq[choose])

plot_flux = flux[choose] * 1e-3 * 1e-23 * freq[choose] * 1e9 / 1e-13
plot_eflux = eflux[choose] * 1e-3 * 1e-23 * freq[choose] * 1e9 / 1e-13
ax.errorbar(
        days[choose], plot_flux, plot_eflux,
        fmt='.', c='k', lw=0.5)
ax.plot(
        days[choose], plot_flux, linestyle='-', c='k',
        label="%s-%s GHz" %(low_freq, max_freq), lw=1.0)

ax.set_ylabel("$\\nu f_{\\nu}$", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)

ax.text(0.05, 0.9, "SMA", transform=ax.transAxes,
        fontsize=14, verticalalignment='top')
ax.text(0.9, 0.9, 
"$10^{-13} \mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-2}$", 
fontsize=14, transform=ax.transAxes, 
verticalalignment='top', horizontalalignment='right')

ax = axarr[1]
ax.errorbar((t-t[0])+3, cts, yerr=ects, fmt='.', c='k')
#ax.plot((t-t[0])+3, cts, c='grey', lw=2.0)
ax.set_xlabel("Days Since 16 June UT", fontsize=14)
ax.set_ylabel(
r"$\nu f_\nu$", fontsize=16)

ax.text(0.9, 0.9, 
        "$10^{-11} \, \mathrm{erg \, s^{-1}} \, \mathrm{cm}^{-2}$",
        fontsize=14, transform=ax.transAxes, 
        verticalalignment='top', horizontalalignment='right')
ax.yaxis.set_tick_params(labelsize=14)
ax.xaxis.set_tick_params(labelsize=14)

#for ii in [12.9, 21.4, 24.4]:
#    ax.axvline(x=ii, alpha=0.5, lw=2)

ax.text(
        0.1, 0.9, "Swift/XRT", 
        transform=ax.transAxes, fontsize=14, verticalalignment='top')

plt.show()
#plt.savefig("xray_vs_radio.png")
