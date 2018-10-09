""" Analysis of the spectral index of AT2018cow """

import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
from lmfit import Model
import emcee
import corner
from get_radio import get_data_all, get_spectrum

def run(nu, f, ef):
    npts = len(nu)
    ndraw = 10000

    # Make 10000 new versions of a new flux array, repopulated with the flux 
    # from the Gaussian dist with STD of that uncertainty

    f_new = np.zeros((npts, ndraw))
    # For each version, do a least squares fit in log space
    # and save the resulting spectral index and y-intercept
    for ii,fval in enumerate(f):
        # log uncertainty is dx/x
        f_new[ii,:] = np.random.normal(loc=fval, scale=ef[ii]/fval, size=ndraw)

    ms = np.zeros(ndraw)
    bs = np.zeros(ndraw)

    for ii in np.arange(ndraw-1):
        m,b = np.polyfit(np.log10(nu), np.log10(f_new[:,ii]), deg=1)
        ms[ii] = m
        bs[ii] = b

    # Measure the mean and standard deviation
    mfit = np.mean(ms)
    emfit = np.std(ms)
    bfit = np.mean(bs)

    # Plot the result
    fig,ax = plt.subplots(1,1)
    plt.errorbar(nu, f, yerr=ef, c='k', fmt='.')
    plt.yscale('log')
    plt.xscale('log')
    xfit = np.logspace(np.log10(min(nu)), np.log10(max(nu)))
    yfit = 10**(mfit*np.log10(xfit)+bfit)
    plt.plot(xfit, yfit, c='r', ls='--')
    plt.text(0.1, 0.9, "m=%s +/- %s" %(mfit, emfit), transform=ax.transAxes)
    plt.text(0.1, 0.8, "b=%s" %(bfit), transform=ax.transAxes)

    plt.show()

    return mfit, emfit, bfit


def run_day_10():
    tel, freq, day, flux, eflux_form, eflux_sys = get_data_all()
    print("Analysis of Day 10" %day)
    choose_tel = tel == 'SMA'
    choose_freq = freq > 100
    choose_day = np.logical_and(day >= 10.2, day <= 10.5)
    choose = choose_tel & choose_freq & choose_day
    nu = freq[choose]
    print("Frequency")
    print(nu)
    f = flux[choose]
    print("Flux")
    print(f)
    print("Unc. on Flux")
    #ef = np.sqrt(eflux_form[choose]**2 + eflux_sys[choose]**2)
    ef = 0.15*f
    print(ef)
    run(nu, f, ef)


def run_day_14():
    tel, freq, day, flux, eflux_form, eflux_sys = get_data_all()
    choose_day = np.logical_and(
            day==14,
            np.logical_or(tel=='SMA', tel=='ALMA'))
    print("Analysis of Day 14")
    choose = np.logical_and(choose_day, freq > 110)
    nu = freq[choose]
    print("Frequency")
    print(nu)
    f = flux[choose]
    print("Flux")
    print(f)
    print("Unc. on Flux")
    ef = 0.1*flux[choose]
    choose_sma = tel[choose]=='SMA'
    ef[choose_sma] = 0.2*f[choose_sma]
    print(ef)
    run(nu, f, ef)


def run_day_22():
    freq,flux = get_spectrum(22)
    print("Analysis of Day 22")
    choose = freq > 110
    nu = freq[choose]
    print("Frequency")
    print(nu)
    f = flux[choose]
    print("Flux")
    print(f)
    print("Unc. on Flux")
    ef = 0.1*flux[choose]
    ef[nu>200] = 0.2*f[nu>200]
    print(ef)
    run(nu, f, ef)


if __name__=="__main__":
    run_day_10()
    #run_day_14()
    #run_day_22()
