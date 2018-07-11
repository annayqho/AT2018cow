""" Make a grid of light curves from the SMA, ATCA, and Swift/XRT """

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
import numpy as np
from astropy.table import Table

dat = Table.read(
    "../data/radio_lc.dat", delimiter="&", format='ascii.no_header')
tel = np.array(dat['col2'])
choose = tel == 'SMA'

days = np.array(dat['col1'][choose])
freq = np.array(dat['col3'][choose]).astype(float)
flux_raw = np.array(dat['col4'][choose])
flux = np.array(
        [float(val.split("pm")[0][1:]) for val in flux_raw])
flux_err = np.array(
        [float(val.split("pm")[1][0:-1]) for val in flux_raw])

choose = np.logical_and(freq > 230, freq < 235)
low_freq = min(freq[choose])
max_freq = max(freq[choose])

plt.errorbar(
        days[choose], flux[choose], yerr=flux_err[choose], 
        mfc='white', mec='black', 
        fmt='.', marker='o', label="%s-%s GHz" %(low_freq, max_freq), c='k')

choose = np.logical_and(freq > 240, freq < 245)
low_freq = min(freq[choose])
max_freq = max(freq[choose])

plt.errorbar(
        days[choose], flux[choose], yerr=flux_err[choose], 
        mfc='white', mec='black', 
        fmt='.', marker='s', label="%s-%s GHz" %(low_freq, max_freq), c='k')

choose = np.logical_and(freq > 230, freq < 245)
plt.plot(
        days[choose], flux[choose], linestyle='--', c='k')

choose = np.logical_and(freq > 330, freq < 334)
low_freq = min(freq[choose])
max_freq = max(freq[choose])

plt.errorbar(
        days[choose], flux[choose], yerr=flux_err[choose], 
        mfc='white', mec='black', ms=8,
        fmt='.', marker='v', label="%s-%s GHz" %(low_freq, max_freq), c='k')

choose = np.logical_and(freq > 357, freq < 361)
low_freq = min(freq[choose])
max_freq = max(freq[choose])

plt.errorbar(
        days[choose], flux[choose], yerr=flux_err[choose], 
        mfc='white', mec='black', ms=8,
        fmt='.', marker='^', label="%s-%s GHz" %(low_freq, max_freq), c='k')

choose = np.logical_or(
        np.logical_and(freq > 330, freq < 334),
        np.logical_and(freq > 357, freq < 361))

plt.plot(
        days[choose], flux[choose], linestyle='-', c='k')

plt.xlabel("Days Since Discovery", fontsize=16)
plt.ylabel("$F_{\\nu}$ [mJy]", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.legend()

plt.tight_layout()
plt.show()
