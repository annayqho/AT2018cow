""" Analysis of the spectral index of AT2018cow """

import numpy as np
from get_radio import get_data_all

tel, freq, day, flux, eflux_form, eflux_sys = get_data_all()

choose = np.logical_and(tel=='SMA', day == 7)

nu = freq[choose]
print(nu)
f = flux[choose]
print(f)
ef = eflux_form[choose]
print(ef)

npts = len(nu)

# if there are fewer than 2 values, don't bother
if npts < 2:
    print("fewer than two points, can't calculate an index")

# if there are two values
if npts == 2:
    alpha = np.log10(f[-1]/f[0])/np.log10(nu[-1]/nu[0])
    ealpha = (1/np.log10(nu[-1]/nu[0])) * (ef[-1]/f[-1]-ef[0]/f[0])
    print(alpha, ealpha)

# if there are more than two values, do a fit to the data
if npts > 2:
    out = np.polyfit(
            np.log10(nu), np.log10(f), w=1/ef, deg=1, full=True, cov=True)
    # according to the np.polyfit documentation,
    # for gaussian uncertainties, use 1/sigma not 1/sigma**2
    print(out)
