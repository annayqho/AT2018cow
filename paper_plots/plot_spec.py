# before running, don't forget to
# delete the lines with AMI
# run :%s///g

import matplotlib
from matplotlib import rc
rc("font", family="serif")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15

if __name__=="__main__":
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

    # plot the ALMA Day 22 and SMA Day 24 data by itself as well
    choose_alma = np.logical_and(days == 22, tel == 'ALMA')
    choose_sma = np.logical_and(days == 24, tel == 'SMA')
    choose = np.logical_or(choose_alma, choose_sma)
    order = np.argsort(freq[choose])

    ax = axarr[1]
    ax.errorbar(
            freq[choose][order], flux[choose][order], 
            yerr=flux_err[choose][order],
            mec='black', fmt='.', c='k')
    # temporary ATCA point
    ax.errorbar(34, 12, yerr=0.1*12, mec='black', fmt='.', c='k')
    ax.text(
            0.9, 0.1, "ATCA Day 19, ALMA Day 22, SMA Day 24", 
            transform=ax.transAxes, horizontalalignment='right', fontsize=14)

    # show a fit
    x = np.linspace(34, 90)
    y = 13*(x/34)**2
    plt.plot(x, y, c='black', lw=1, ls='-')

    # shallower part
    x = np.linspace(90, 120)
    y = 96*(x/103)**(1/3)
    plt.plot(x, y, c='black', lw=1, ls='-')

    # power-law index
    x = np.linspace(120, 350)
    alpha = -1.11
    y = 80*((x/150)**alpha)
    plt.plot(x, y, c='black', lw=1, ls='-')

    # show critical frequencies
    plt.axvline(x=90, c='grey', ls='--', alpha=0.8)
    plt.axvline(x=120, c='grey', ls='--', alpha=0.8)

    ax.tick_params(axis='both', labelsize=14)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlabel("Frequency [GHz]", fontsize=16)
    ax.set_ylabel("Flux [mJy]", fontsize=16)
    ax.set_xscale('log')
    ax.set_yscale('log')

    #plt.show()
    plt.savefig("spec.png")
