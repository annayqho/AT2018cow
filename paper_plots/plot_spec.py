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

def spec_sequence(days, freq, flux, flux_err, lims):
    colors = ['pink', 'red', 'orange', 'yellow', 'green', 'blue', 'purple', 'black', 'pink', 'red', 'orange', 'yellow', 'green', 'blue', 'purple', 'black']

    for ii,day in enumerate(days):
        # plot detections
        print(day)
        choose = np.where(np.logical_and(t==day, lims==False))[0]
        #plt.errorbar(
        #        freq[choose], flux[choose], yerr=flux_err[choose], 
        #        color=colors[ii], fmt='.')
        if ii < 8:
            plt.plot(
                    freq[choose], flux[choose], color=colors[ii], 
                    lw=2.0, ls='--')
        else:
            plt.plot(
                    freq[choose], flux[choose], color=colors[ii], 
                    lw=2.0, ls='-')
        plt.show()


        # plot non-detections
        #choose = np.where(np.logical_and(t==day, lims==True))[0]
        #plt.scatter(
        #        freq[choose], flux[choose]+1, color=colors[ii], marker='v')

    # formatting
    #plt.xlim(0.1, 360)
    #plt.ylim(0.1, 42)
    plt.ylabel("Flux [mJy]", fontsize=12)
    plt.xscale('log')
    plt.yscale('log')
    # plt.text(
    #         0.10, 0.90, "$Delta t$ = %s d" %day, 
    #         horizontalalignment='left', verticalalignment='top',
    #         transform=axarr[ii].transAxes)

    plt.xlabel("Freq [GHz]", fontsize=12)
    plt.tight_layout()
    #plt.show()
    #plt.savefig("spec_sequence.png")


if __name__=="__main__":
    dat = Table.read(
        "../data/radio_lc.dat", delimiter="&", format='ascii.no_header')
    tel = np.array(dat['col2'])
    choose = np.array([val in np.array(['SMA', 'ATCA', 'ALMA']) for val in tel])

    tel = tel[choose]
    days = np.array(dat['col1'][choose])
    freq = np.array(dat['col3'][choose]).astype(float)
    flux_raw = np.array(dat['col4'][choose])
    flux = np.array(
            [float(val.split("pm")[0][1:]) for val in flux_raw])
    flux_err = np.array(
            [float(val.split("pm")[1][0:-1]) for val in flux_raw]) 

    # first, try to show the spectrum on one plot
    fig = plt.figure(figsize=(8,6))

    # plot the spectrum from Day 10 data
    choose = days == 10
    order = np.argsort(freq[choose])
    plt.errorbar(
            freq[choose][order], flux[choose][order], 
            yerr=flux_err[choose][order],
            mfc='white', mec='black', fmt=':', marker='o',
            label="$\Delta t = 10$", c='k')

    # plot the spectrum from Day 13-14 data
    # ATCA from Day 13, SMA and ALMA from Day 14
    choose = np.logical_or(
            np.logical_and(days == 13, tel == 'ATCA'),
            np.logical_and(days == 14, tel == 'SMA'))

    order = np.argsort(freq[choose])
    plt.errorbar(
            freq[choose][order], flux[choose][order], 
            yerr=flux_err[choose][order],
            mfc='white', mec='black', fmt='--', marker='s',
            label="$\Delta t = 13-14$", c='k')

    # plot the spectrum from the Day 22 data
    choose = np.logical_and(days >= 22, days <= 24)

    order = np.argsort(freq[choose])
    plt.errorbar(
            freq[choose][order], flux[choose][order], 
            yerr=flux_err[choose][order],
            mfc='white', mec='black', fmt='-', marker='v', ms=8, ls='-',
            label="$\Delta t = 22-24$", c='k')

    #spec_sequence(day_bins, freq, flux, flux_err, lims)

    #spectral_index(t, freq, flux, flux_err, lims)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel("Frequency [GHz]", fontsize=16)
    plt.ylabel("Flux [mJy]", fontsize=16)

    plt.legend()

    plt.show()
