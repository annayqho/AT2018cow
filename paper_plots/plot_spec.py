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


def synchrotron(x, nu_a, F_a, nu_m, F_m, alpha):
    """ Fitting function """
    y = np.zeros(len(x))
    choose = x <= nu_a
    y[choose] = F_a * (x[choose]/nu_a)**2
    choose = np.logical_and(x > nu_a, x <= nu_m)
    y[choose] = F_a * (x[choose]/nu_a)**(1/3)
    choose = x > nu_m
    y[choose] = F_m * (x[choose]/nu_m)**alpha
    return y


def fit_spec(freq, flux, flux_err):
    """ Fit the spectrum on Day 22 to a synchrotron spectrum """
    nu_a = 90
    F_a = 90
    nu_m = 100
    F_m = 90
    alpha = -1.0
    p0 = [nu_a, F_a, nu_m, F_m, alpha]

    popt, pcov = curve_fit(
            synchrotron, freq, flux, p0 = p0,
            sigma = flux_err, absolute_sigma=True)

    xfit = np.linspace(min(freq), max(freq))
    yfit = synchrotron(xfit, *popt)




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


if __name__=="__main__":
    tel, freq, days, flux, flux_err = get_data()

    fig, axarr = plt.subplots(2, 1, figsize=(8,8))

    # top panel: evolution of spectra
    ax = axarr[0]

    # plot the spectrum from Day 10 data
    choose = days == 10
    order = np.argsort(freq[choose])
    det = flux_err[choose][order] > 0
    ax.errorbar(
            freq[choose][order][det], flux[choose][order][det], 
            yerr=flux_err[choose][order][det],
            mfc='white', mec='black', fmt=':', marker='o',
            label="$\Delta t = 10$", c='k')

    # plot the spectrum from Day 13-14 data
    # ATCA from Day 13, SMA and ALMA from Day 14
    choose = np.logical_or(
            np.logical_and(days == 13, tel == 'ATCA'),
            np.logical_and(days == 14, tel != 'ATCA'))

    order = np.argsort(freq[choose])

    # for detections
    det = flux_err[choose][order] > 0
    ax.errorbar(
            freq[choose][order][det], flux[choose][order][det], 
            yerr=flux_err[choose][order][det],
            mfc='white', mec='black', fmt='--', marker='s',
            label="$\Delta t = 13, 14$", c='k')

    # plot the spectrum from the Day 22 data
    choose = days == 22
    choose = np.logical_or(days == 22, days == 24)

    order = np.argsort(freq[choose])
    ax.errorbar(
            freq[choose][order], flux[choose][order], 
            yerr=flux_err[choose][order],
            mfc='white', mec='black', fmt='-', marker='v', ms=8, ls='-',
            label="$\Delta t = 22, 24$", c='k')

    ax.tick_params(axis='both', labelsize=14)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_ylabel("Flux [mJy]", fontsize=16)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.legend()


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

    ax = axarr[1]
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

    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlabel("Frequency [GHz]", fontsize=16)
    ax.set_ylabel("Flux [mJy]", fontsize=16)
    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.show()
    #plt.savefig("spec.png")
