""" Analysis of the spectral index of AT2018cow """

import numpy as np
from scipy.optimize import curve_fit
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
    xfit = np.linspace(min(nu), max(nu), 100)
    yfit = f[0] * 10**(alpha*np.log10(xfit/nu[0]))
    plt.scatter(nu, f, c='k')
    plt.plot(xfit, yfit, ls=':', c='r')

# if there are more than two values, do a fit to the data
if npts > 2:
    logx = np.log10(nu)
    logy = np.log10(f)
    logyerr = ef/y
    powerlaw = lambda x, amp, index: amp * (x**index)
    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
    pinit = [1.0, -1.0]
    out = optimize.leastsq(errfunc, pinit,
            args=(logx, logy, logyerr), full_output=1)
    pfinal = out[0]
    covar = out[1]
    index = pfinal[1]
    amp = 10**pfinal[0]
    indexErr = np.sqrt( covar[1][1] )
    ampErr = np.sqrt( covar[0][0] ) * amp

    plt.plot(xdata, powerlaw(xdata, amp, index))     # Fit
    plt.errorbar(xdata, ydata, yerr=yerr, fmt='k.')  # Data
     
plt.show()
