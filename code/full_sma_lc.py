""" Full SMA light curve breaking data into all frequency tunings, 
for the Appendix """

import numpy as np
import matplotlib
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/code")
from scale_fluxes import sma_lc
from get_radio import get_data_all

out = sma_lc()
t, nu, f, ef = out[0]

# A non-detection will either have a negative flux value
# or an uncertainty that is larger in magnitude than the flux value
det = np.logical_and(f > 0, np.abs(f) > ef)

# Only plot detections
t = t[det]
nu = nu[det]
f = f[det]
ef = ef[det]

npts = len(np.unique(nu))

c = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c']
markers = ['o', 's', 'v', '^']
nlab = len(c)

fig,axarr = plt.subplots(4,1, figsize=(8,12), sharex=True, sharey=True)
plt.subplots_adjust(hspace=0.1)

ax = axarr[0]
ax.text(
        0.05, 0.95, "\\textbf{193.5-218 GHz}", fontsize=14, transform=ax.transAxes,
        horizontalalignment='left', verticalalignment='top')
choose = nu < 220

for ii,n in enumerate(np.unique(nu[choose])):
    unc = np.sqrt(ef[nu==n]**2 + (0.15*f[nu==n])**2)
    ax.errorbar(
            t[nu==n], f[nu==n], unc, 
            c=c[ii%nlab], 
            fmt=markers[(int(np.floor(ii/nlab))+ii)%nlab],
            label="%s GHz" %n, ms=7)
    ax.plot(
            t[nu==n], f[nu==n], c=c[ii%nlab], lw=1.0)

ax.legend(ncol=3)
ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)

ax = axarr[1]
ax.text(
        0.05, 0.95, "\\textbf{225-234.6 GHz}", fontsize=14, transform=ax.transAxes,
        horizontalalignment='left', verticalalignment='top')
choose = np.logical_and(nu >= 220, nu < 240)
for ii,n in enumerate(np.unique(nu[choose])):
    unc = np.sqrt(ef[nu==n]**2 + (0.15*f[nu==n])**2)
    ax.errorbar(
            t[nu==n], f[nu==n], unc, 
            c=c[ii%nlab],
            fmt=markers[(int(np.floor(ii/nlab))+ii)%nlab],
            label="%s GHz" %n, ms=7)
    ax.plot(t[nu==n], f[nu==n], c=c[ii%len(c)], lw=1.0)
ax.legend(ncol=3)
ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)

ax = axarr[2]
ax.text(
        0.05, 0.95, "\\textbf{241-283 GHz}", fontsize=14, transform=ax.transAxes,
        horizontalalignment='left', verticalalignment='top')
choose = np.logical_and(nu >= 240, nu < 300)
for ii,n in enumerate(np.unique(nu[choose])):
    unc = np.sqrt(ef[nu==n]**2 + (0.15*f[nu==n])**2)
    ax.errorbar(
            t[nu==n], f[nu==n], unc, 
            c=c[ii%nlab], 
            fmt=markers[(int(np.floor(ii/nlab))+ii)%nlab], 
            label="%s GHz" %n, ms=7)
    ax.plot(t[nu==n], f[nu==n], c=c[ii%len(c)], lw=1.0)
ax.legend(ncol=3)
ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)

ax = axarr[3]
ax.text(
        0.05, 0.95, "\\textbf{330.8-360.8 GHz}", fontsize=14, 
        transform=ax.transAxes,
        horizontalalignment='left', verticalalignment='top')
choose = nu > 300
for ii,n in enumerate(np.unique(nu[choose])):
    unc = np.sqrt(ef[nu==n]**2 + (0.15*f[nu==n])**2)
    ax.errorbar(
            t[nu==n], f[nu==n], unc, 
            c=c[0:3][ii%3],
            fmt=markers[(int(np.floor(ii/nlab))+ii)%nlab],
            label="%s GHz" %n, ms=7)
    ax.plot(t[nu==n], f[nu==n], c=c[ii%len(c)], lw=1.0)
ax.legend(ncol=3, loc='lower right')
ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)

ax.xaxis.set_tick_params(labelsize=14)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1.0,100)
ax.set_xlabel("Time Since June 16 UT (MJD 58285) [d]", fontsize=16)

#plt.show()
plt.savefig("full_lc.png")
