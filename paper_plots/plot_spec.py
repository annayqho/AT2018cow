# before running, don't forget to
# delete the lines with AMI
# run :%s///g

import matplotlib
from matplotlib import rc
rc("font", family="serif")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from astropy.table import Table
from astropy.cosmology import Planck15


def synchrotron(x, nu_a, F_a, alpha):
    """ Fitting function """
    y = np.zeros(len(x))
    nu_m = 100
    F_m = 95
    choose = x <= nu_a
    y[choose] = F_a * (x[choose]/nu_a)**2
    choose = np.logical_and(x > nu_a, x <= nu_m)
    y[choose] = F_a * (x[choose]/nu_a)**(1/3)
    choose = x > nu_m
    y[choose] = F_m * (x[choose]/nu_m)**alpha
    return y


def self_abs(x, b):
    """ Fitting function for self-absorbed part """
    y = 10**(2*np.log10(x)+b)
    return y


def fit_spec(freq, flux, flux_err):
    """ Fit the spectrum on Day 22 to a synchrotron spectrum """
    nu_a = 90
    F_a = 90
    alpha = -1.0
    p0 = [nu_a, F_a, alpha]

    popt, pcov = curve_fit(
            synchrotron, freq, flux, p0 = p0,
            sigma = flux_err, absolute_sigma=True)

    xfit = np.linspace(min(freq), max(freq))
    yfit = synchrotron(xfit, *popt)


def fit_self_abs(freq, flux, flux_err=None):
    """ Fit the self-absorbed part of the spectrum """
    if flux_err:
        popt, pcov = curve_fit(
                self_abs, freq, flux,
                sigma = flux_err, absolute_sigma = True)
    else:
        popt, pcov = curve_fit(
                self_abs, freq, flux)
    return popt[0]



def get_data():
    """ Get the light curves """
    dat = Table.read(
        "../data/radio_lc.dat", delimiter="&", format='ascii.no_header')
    days = np.array(dat['col1'])
    tel = np.array(dat['col2'])
    freq = np.array(dat['col3']).astype(float)
    flux_raw = np.array(dat['col4'])
    flux = np.zeros(len(flux_raw))
    flux_err = np.zeros(len(flux_raw))

    # Now, go through and apply the appropriate uncertainties
    for ii,val in enumerate(flux_raw):
        if tel[ii] == 'ALMA':
            flux[ii] = val.split("$")[1]
            if freq[ii] < 200: # bands 3 and 4
                flux_err[ii] = 0.05*flux[ii]
            elif freq[ii] > 600: # band 9
                flux_err[ii] = 0.20*flux[ii]
            else: # bands 7 and 8
                flux_err[ii] = 0.10*flux[ii]
        elif tel[ii] == 'SMA':
            flux[ii] = float(val.split("pm")[0][1:]) 
            flux_err_temp = float(val.split("pm")[1][0:-1])
            # add 10% uncertainty
            flux_err[ii] = np.sqrt(
                    (0.10*flux[ii])**2 + flux_err_temp**2)
        elif tel[ii] == 'ATCA':
            flux[ii] = float(val.split("pm")[0][1:]) 
            # add 10% uncertainty
            flux_err_temp = float(val.split("pm")[1][0:-1])
            if flux_err_temp < 0:
                flux_err[ii] = -99
            else:
                flux_err[ii] = np.sqrt(
                        (0.10*flux[ii])**2 + flux_err_temp**2)
    return tel, freq, days, flux, flux_err


def run_day(day):
    """ Run this for one day """
    choose = days == day
    order = np.argsort(freq[choose])
    det = flux_err[choose][order] > 0
    ax.errorbar(
            freq[choose][order][det], flux[choose][order][det], 
            yerr=flux_err[choose][order][det], fmt='o',
            label="$\Delta t = %s$" %day, c='k')
    ax.scatter(
            5.5, 0.15, marker='v', c='k')

    b = fit_self_abs(
            freq[choose][order][det][0:2],
            flux[choose][order][det][0:2],
            flux_err[choose][order][det][0:2])
    xfit = np.linspace(5, 200)
    yfit = 10**(2*np.log10(xfit)+b)
    ax.plot(xfit, yfit, ls=':', c='k')

    ax.text(
            0.4, 0.45, "$F_\\nu \propto \\nu^2$", 
            transform=ax.transAxes, horizontalalignment='right', fontsize=14)

    # fit and plot a nu^something line
    m,b = np.polyfit(
            np.log10(freq[choose][order][det][2:]),
            np.log10(flux[choose][order][det][2:]),
            deg = 1,
            w=1/flux_err[choose][order][det][2:]**2)
    print(m)
    xfit = np.linspace(80, 500)
    yfit = 10**(m*np.log10(xfit)+b)
    ax.plot(xfit, yfit, ls=':', c='k')
    ax.text(
            0.8, 0.85, "$F_\\nu \propto \\nu^{-1.8}$", 
            transform=ax.transAxes, horizontalalignment='left', fontsize=14)

    # the lines join at (144, 92)
    ax.text(
            0, 0.9, "$(\\nu_p, F_p) = (144\,\mathrm{GHz}, 92\,\mathrm{mJy}$)", 
            transform=ax.transAxes, horizontalalignment='left', fontsize=14)
    ax.text(
            0.9, 0.1, "ATCA & SMA, Day %s" %day, 
            transform=ax.transAxes, horizontalalignment='right', fontsize=14)

    ax.tick_params(axis='both', labelsize=14)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_ylabel("Flux [mJy]", fontsize=16)
    ax.set_xscale('log')
    ax.set_yscale('log')


if __name__=="__main__":
    tel, freq, days, flux, flux_err = get_data()

    fig, axarr = plt.subplots(
            3, 1, figsize=(8,8), sharex=True, sharey=True)

    # top panel: evolution of spectra
    ax = axarr[0]

    # ATCA from Day 13, SMA and ALMA from Day 14
    ax = axarr[1]

    choose = np.logical_or(
            np.logical_and(days == 13, tel == 'ATCA'),
            np.logical_and(days == 14, tel != 'ATCA'))

    order = np.argsort(freq[choose])

    # for detections
    det = flux_err[choose][order] > 0
    ax.errorbar(
            freq[choose][order][det], flux[choose][order][det], 
            yerr=flux_err[choose][order][det],
            label="$\Delta t = 13, 14$", fmt='o', c='k')

    # plot a nu^2 line
    b = fit_self_abs(
            freq[choose][order][det][0:5],
            flux[choose][order][det][0:5],
            flux_err[choose][order][det][0:5])
    xfit = np.linspace(5, 150)
    yfit = 10**(2*np.log10(xfit)+b)
    ax.plot(xfit, yfit, ls=':', c='k')
    ax.text(
            0.4, 0.5, "$F_\\nu \propto \\nu^2$", 
            transform=ax.transAxes, horizontalalignment='right', fontsize=14)

    # fit and plot a nu^something line
    m,b = np.polyfit(
            np.log10(freq[choose][order][det][5:]),
            np.log10(flux[choose][order][det][5:]),
            deg = 1,
            w=1/flux_err[choose][order][det][5:]**2)
    print(m)
    xfit = np.linspace(80, 500)
    yfit = 10**(m*np.log10(xfit)+b)
    ax.plot(xfit, yfit, ls=':', c='k')
    ax.text(
            0.9, 0.85, "$F_\\nu \propto \\nu^{-0.8}$", 
            transform=ax.transAxes, horizontalalignment='right', fontsize=14)

    # Explain which day and instruments this is
    ax.text(
            0, 0.9, "$(\\nu_p, F_p) = (111\,\mathrm{GHz}, 73\,\mathrm{mJy}$)", 
            transform=ax.transAxes, horizontalalignment='left', fontsize=14)
    ax.text(
            0.9, 0.1, "ATCA on Day 13, ALMA & SMA on Day 14", 
            transform=ax.transAxes, horizontalalignment='right', fontsize=14)

    ax.tick_params(axis='both', labelsize=14)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_ylabel("Flux [mJy]", fontsize=16)
    ax.set_xscale('log')
    ax.set_yscale('log')

    # bottom panel: a fit to the spectrum on Day 22
    # ALMA data
    choose_alma = np.logical_and(days == 22, tel == 'ALMA')

    # interpolate the SMA data at 215 and 230 GHz
    choose_days = np.logical_or(days==20, days==24)
    choose_temp = np.logical_and(tel == 'SMA', freq == 215.5)
    choose_215 = np.logical_and(choose_days, choose_temp)
    choose_temp = np.logical_and(tel == 'SMA', freq == 231.5)
    choose_231 = np.logical_and(choose_days, choose_temp)

    # take the average
    f_215 = np.average(
            flux[choose_215], weights=1/flux_err[choose_215]**2)
    f_231 = np.average(
            flux[choose_231], weights=1/flux_err[choose_231]**2)

    #choose = np.logical_or(choose_alma, choose_sma)
    choose = choose_alma
    order = np.argsort(freq[choose])

    ax = axarr[2]
    ax.errorbar(
            freq[choose][order], flux[choose][order], 
            yerr=flux_err[choose][order],
            mec='black', fmt='o', c='k')
    # temporary ATCA point
    ax.errorbar(34, 12, yerr=0.2*12, mec='black', fmt='o', c='k')
    # temporary SMA points
    print(f_215, f_231)
    ax.errorbar(
            [215.5, 231.5], [f_215, f_231], yerr=0.1*np.array([f_215, f_231]), 
            mec='black', fmt='o', c='k')
    ax.text(
            0.9, 0.1, "ATCA, ALMA, & SMA, Day 22", 
            transform=ax.transAxes, horizontalalignment='right', fontsize=14)
    ax.text(
            0, 0.9, "$(\\nu_p, F_p) = (92\,\mathrm{GHz}, 92\,\mathrm{mJy}$)", 
            transform=ax.transAxes, horizontalalignment='left', fontsize=14)

    ax.text(
            0.6, 0.8, "$F_\\nu \propto \\nu^2$", 
            transform=ax.transAxes, horizontalalignment='right', fontsize=14)

    # fit a nu^2 line to the ATCA point and the lowest-freq ALMA point
    x = np.array([34, freq[choose][order][0]])
    y = np.array([12, flux[choose][order][0]])
    yerr = np.array([0.1*12, flux_err[choose][order][0]])
    b = fit_self_abs(x, y, yerr)
    xfit = np.linspace(20, 150)
    yfit = 10**(2*np.log10(xfit)+b)
    plt.plot(xfit, yfit, c='k', ls=':')

    # fit and plot a nu^something line
    freq_fit = np.hstack((freq[choose][order], 215, 231))
    flux_fit = np.hstack((flux[choose][order], f_215, f_231))
    flux_err_fit = np.hstack((flux_err[choose][order], 0.1*f_215, 0.1*f_231))
    m,b = np.polyfit(
            np.log10(freq_fit[4:]), np.log10(flux_fit[4:]), deg = 1,
            w=1/flux_err_fit[4:]**2)
    print(m)
    xfit = np.linspace(80, 500)
    yfit = 10**(m*np.log10(xfit)+b)
    ax.plot(xfit, yfit, ls=':', c='k')
    ax.text(
            0.95, 0.85, "$F_\\nu \propto \\nu^{-1.2}$", 
            transform=ax.transAxes, horizontalalignment='right', fontsize=14)

    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlabel("Frequency [GHz]", fontsize=16)
    ax.set_ylabel("Flux [mJy]", fontsize=16)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(0.1,300)

    #plt.show()
    plt.savefig("spec.png")
