import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.cosmology import Planck15
import sys
sys.path.append('/Users/annaho/Dropbox/Projects/Research/AT2018cow/code')
from scale_fluxes import sma_lc

a, b = sma_lc()
dt, f, ef = b
ef_comb = np.sqrt(ef**2 + (0.15*f)**2)
plt.scatter(
        dt, f,
        marker='s', c='k',
        label="230.6--234.6 GHz")
plt.errorbar(
        dt, f, ef_comb, 
        fmt='s', c='k', lw=1.5)
plt.plot(
        dt, f, linestyle='-', c='k', lw=1.5)
plt.ylabel("$f_{\\nu}$ [mJy]", fontsize=16)
plt.xlabel("$\Delta t$ from 16 June UT [days]", fontsize=16)
plt.tick_params(axis='both', labelsize=14)
plt.legend(fontsize=12, loc='lower left')
plt.ylim(3,70)
plt.yscale('log')
plt.xscale('log')
plt.gca().get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.gca().get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.yticks([5,10,20,50])
plt.xticks([5,10,20,60])
plt.show()
#plt.savefig("lc_230ghz.png")
