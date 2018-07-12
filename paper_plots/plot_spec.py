# before running, don't forget to
# delete the lines with AMI
# run :%s/\\//g

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
    #         0.10, 0.90, "$\Delta t$ = %s d" %day, 
    #         horizontalalignment='left', verticalalignment='top',
    #         transform=axarr[ii].transAxes)

    plt.xlabel("Freq [GHz]", fontsize=12)
    plt.tight_layout()
    #plt.show()
    #plt.savefig("spec_sequence.png")


if __name__=="__main__":
    dat = Table.read("../data/radio_lc.dat", delimiter="&", format='ascii')

    t = dat['dt']
    freq = dat['Frequency']
    flux_raw = dat['Flux']

    # put in order of freq
    order = np.argsort(freq)
    t = t[order]
    freq = freq[order]
    flux_raw = flux_raw[order]

    npts = len(flux_raw)
    flux = np.zeros(npts)
    flux_err = np.zeros(npts)
    lims = np.zeros(npts, dtype=bool)

    # formatting
    for ii,val in enumerate(flux_raw):
        if "<" in val:
            flux[ii] = val.split("$")[1][1:]
            lims[ii] = True
        else:
            temp = val.split("pm")
            flux[ii] = float(temp[0][1:])
            flux_err[ii] = float(temp[1][0:-1])
            lims[ii] = False

    # unique days
    day_bins = np.unique(t)
    # cols = ['red', 'orange', 'yellow', 'green', 'blue', 'purple', 'black']
    # cols = ['#f2f0f7', '#dadaeb', '#bcbddc', '#9e9ac8', '#807dba', '#6a51a3', '#4a1486']
    #spec_sequence(day_bins, freq, flux, flux_err, lims)
    light_curve(t, freq, flux, flux_err)

    #spectral_index(t, freq, flux, flux_err, lims)
