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


def light_curve(t, freq, flux, flux_err):
    #conv = 1e-3 * 10**(-23)
    #d = Planck15.luminosity_distance(z=0.014).cgs.value
    #lum = flux * conv * 4 * np.pi * d**2
    #lum_err = flux_err * conv * 4 * np.pi * d**2
    freq_bins = [5.5, 9, 34, 230, 340]
    window = [1, 2, 2, 20, 20]
    marker_type = ['*', 'o', 's', '^', 'v']
    line_style = ['--', '-', ':', '--', '-.']

    for ii, center_freq in enumerate(freq_bins):
        choose = np.where(np.abs(freq-center_freq) <= window[ii])[0]

        # choose one instance of every day
        select,ind = np.unique(t[choose], return_index=True)

        plt.errorbar(
                t[choose][ind], flux[choose][ind], 
                yerr=flux_err[choose][ind],
                mfc='white', mec='black', fmt='.', marker=marker_type[ii], 
                linestyle=line_style[ii], label="%s GHz" %center_freq, c='k')

    # put the 1998bw point on there
    #plt.scatter(
    #        12, 6.7e28/1e28, marker='X', label="1998bw, SCUBA $150\,$GHz", c='k')

    # make it pretty
    plt.xlabel("Days Since Discovery", fontsize=16)
    plt.ylabel("$F_{\\nu}$ [mJy]", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.yscale('lin')
    
    plt.legend()
    plt.tight_layout()

    #plt.show()
    plt.savefig("light_curve.png")


def plot_spectral_index(d, f, lc, lc_err, is_lim, ind_1, ind_2, mtype, lstyle):
    # Find the overlap in days
    d_choose = np.intersect1d(d[ind_1], d[ind_2])
    choose_1 = np.array(
            [np.where(d[ind_1]==d_val)[0][0] for d_val in d_choose])
    lc_1 = lc[ind_1][choose_1]
    lim_1 = is_lim[ind_1][choose_1]
    choose_2 = np.array(
            [np.where(d[ind_2]==d_val)[0][0] for d_val in d_choose])
    lc_2 = lc[ind_2][choose_2]
    lim_2 = is_lim[ind_2][choose_2]

    # Calculate the spectral index at each epoch
    alpha = np.log10( lc_1 / lc_2 ) / np.log10( f[ind_2] / f[ind_1] )

    # Plot
    plt.errorbar(
            d_choose, alpha,
            mfc='white', mec='black', fmt='.', marker=mtype, 
            linestyle=lstyle, c='k',
            label="%s to %s GHz" %(f[ind_1], f[ind_2]))
    plt.plot(
            d_choose, alpha, color='k', lw=1.0)


def spectral_index(t, freq, flux, flux_err, lims):
    freq_bins = [9, 34, 230, 340]
    d = []
    lc = []
    lc_err = []
    is_lim = []

    # some stuff for plotting
    marker_type = ['o', 's', 'v']
    line_style = ['-', ':', '--']

    for ii, center_freq in enumerate(freq_bins):
        choose = np.where(np.abs(freq-center_freq) <= 0.2*center_freq)[0]

        # choose one instance of every day
        select,ind = np.unique(t[choose], return_index=True)
        d.append(np.array(t[choose][ind]))
        lc.append(np.array(flux[choose][ind]))
        lc_err.append(np.array(flux_err[choose][ind]))
        is_lim.append(np.array(lims[choose][ind]))

    lc = np.array(lc)
    lc_err = np.array(lc_err)
    is_lim = np.array(is_lim)

    plot_spectral_index(
            d, freq_bins, lc, lc_err, is_lim, 
            0, 1, mtype = marker_type[0], lstyle = line_style[0])

    plot_spectral_index(
            d, freq_bins, lc, lc_err, is_lim, 
            1, 2, mtype = marker_type[1], lstyle = line_style[1])

    plot_spectral_index(
            d, freq_bins, lc, lc_err, is_lim, 
            2, 3, mtype = marker_type[2], lstyle = line_style[2])



    plt.xlabel("Days Since Discovery", fontsize=16)
    plt.ylabel(
            "Spectral Index $\\alpha$, $F_\\nu \propto \\nu^{-\\alpha}$", 
            fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
    plt.legend()
    plt.tight_layout()

    #plt.show()
    plt.savefig("spectral_index.png")


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
