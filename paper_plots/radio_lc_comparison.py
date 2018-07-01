""" Compare the radio light curve at ATCA frequencies (say, ATCA at 5.5 
and 9 GHz to those of other potentially comparable events

For now let's just do it at 9 GHz
"""

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
import numpy as np
from astropy.cosmology import Planck15
from astropy.table import Table


fig,axarr = fig.subplots(2,1)

# for AT2018cow
d = Planck15.luminosity_distance(z=0.014).cgs.value
lum = (0.46 * 1e-3 * 10**(-23)) * 4 * np.pi * d**2
lum_err = (0.06 * 1e-3 * 10**(-23)) * 4 * np.pi * d**2
axarr[0].errorbar(10, lum, yerr=lum_err, marker='o', c='k', label="AT2018cow")


# for 1998bw
dat = Table.read(
    "../data/radio_compilations/1998bw.dat", 
    delimiter="&", format='ascii.fixed_width')
freq = dat['freq']
choose = freq==8.64
t = dat['dt'][choose]
flux = dat['flux'][choose] * 1e-3 * 10**(-23)
d = Planck15.luminosity_distance(z=0.0085).cgs.value 
lum = flux * 4 * np.pi * d**2
axarr[0].plot(t, lum, c='grey', label="1998bw")


# for 2009bb
dat = Table.read(
    "../data/radio_compilations/2009bb.dat", 
    delimiter="&", format='ascii.fixed_width')
freq = dat['freq']
choose = freq==8.46
t = dat['dt'][choose]
flux = dat['flux'][choose] * 1e-3 * 10**(-23) 
d = Planck15.luminosity_distance(z=0.0085).cgs.value 
lum = flux * 4 * np.pi * d**2 
axarr[0].plot(t, lum, c='grey', label="2009bb", linestyle='--')

axarr[0].set_ylabel("$L_{9\,\mathrm{GHz}}$ [erg/cm$^2$/s/Hz]", fontsize=16)
axarr[0].set_yscale('log')


### Now again, except this time for mm observations...


# for 1998bw
Planck15.luminosity_distance(z=0.0085).cgs.value 
lum = 39 * 1e-3 * 10**(-23)  * 4 * np.pi * d**2
plt.scatter(12, lum, marker='s', label="1998bw, SCUBA, 150\,GHz")

d = Planck15.luminosity_distance(z=0.014).cgs.value 
lum = (0.46 * 1e-3 * 10**(-23)) * 4 * np.pi * d**2 
lum_err = (0.06 * 1e-3 * 10**(-23)) * 4 * np.pi * d**2 
axarr[0].errorbar(10, lum, yerr=lum_err, marker='o', c='k', label="AT2018cow")

plt.xlabel("Days Since Explosion", fontsize=16)
plt.legend()
plt.tight_layout()

plt.show()
