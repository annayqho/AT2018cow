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
npts = len(np.unique(nu))

c = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#b15928', 'k']

fig,axarr = plt.subplots(4,1, figsize=(8,12), sharex=True, sharey=True)
plt.subplots_adjust(hspace=0)

ax = axarr[0]
choose = nu < 220
for ii,n in enumerate(np.unique(nu[choose])):
    unc = np.sqrt(ef[nu==n]**2 + (0.15*f[nu==n])**2)
    ax.errorbar(t[nu==n], f[nu==n], unc, c=c[ii%len(c)], fmt='.')
    ax.plot(t[nu==n], f[nu==n], c=c[ii%len(c)], label="%s GHz" %n)
ax.legend(ncol=3)
ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)

ax = axarr[1]
choose = np.logical_and(nu >= 220, nu < 240)
for ii,n in enumerate(np.unique(nu[choose])):
    unc = np.sqrt(ef[nu==n]**2 + (0.15*f[nu==n])**2)
    ax.errorbar(t[nu==n], f[nu==n], unc, c=c[ii%len(c)], fmt='.')
    ax.plot(t[nu==n], f[nu==n], c=c[ii%len(c)], label="%s GHz" %n)
ax.legend(ncol=3)
ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)

ax = axarr[2]
choose = np.logical_and(nu >= 240, nu < 300)
for ii,n in enumerate(np.unique(nu[choose])):
    unc = np.sqrt(ef[nu==n]**2 + (0.15*f[nu==n])**2)
    ax.errorbar(t[nu==n], f[nu==n], unc, c=c[ii%len(c)], fmt='.')
    ax.plot(t[nu==n], f[nu==n], c=c[ii%len(c)], label="%s GHz" %n)
ax.legend(ncol=3)
ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)

ax = axarr[3]
choose = nu > 300
for ii,n in enumerate(np.unique(nu[choose])):
    unc = np.sqrt(ef[nu==n]**2 + (0.15*f[nu==n])**2)
    ax.errorbar(t[nu==n], f[nu==n], unc, c=c[ii%len(c)], fmt='.')
    ax.plot(t[nu==n], f[nu==n], c=c[ii%len(c)], label="%s GHz" %n)
ax.legend(ncol=3, loc='lower right')
ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)
ax.yaxis.set_tick_params(labelsize=14)


ax.xaxis.set_tick_params(labelsize=14)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(2,100)
ax.set_xlabel("Time Since June 16 UT (MJD 58285) [d]", fontsize=16)

#plt.show()
plt.savefig("full_lc.png")
