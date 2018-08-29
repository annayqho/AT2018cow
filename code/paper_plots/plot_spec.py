import matplotlib
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/code")
from get_radio import get_data_all, get_spectrum
from synchrotron_fit import self_abs, fit_self_abs


def plot_day(ax,day,nu,flux,islim,formatting=True,fit_peak=True,quad=False):
    """ Run this for one day """
    ax.scatter(
            nu[~islim], flux[~islim], 
            marker='o',label="$\Delta t = %s$" %day, c='k')
    ax.scatter(
            nu[islim], flux[islim], marker='v', c='k')

    if formatting:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylabel("Flux [mJy]", fontsize=16)
        ax.text(
                0.3, 0.45, "$F_\\nu \propto \\nu^2$", 
                transform=ax.transAxes, horizontalalignment='right', 
                fontsize=14)
    if fit_peak:
        choose = np.logical_and(nu<70, ~islim)
        b = fit_self_abs(nu[choose],flux[choose],0.1*flux[choose])
        xfit = np.linspace(5, 200)
        yfit = 10**(2*np.log10(xfit)+b)
        ax.plot(xfit, yfit, ls=':', c='k')

        # fit and plot a nu^something line
        xfit = np.linspace(50, 1000)
        if day == 10:
            m = -1.87
            em = 0.46
            b = 6.0087
        elif day == 14:
            m = -0.93
            em = 0.27
            b = 3.836
        elif day == 22:
            m = -1.1
            em = 0.2
            b = 4.293
        yfit = 10**(m*np.log10(xfit)+b)
        ax.plot(xfit, yfit, ls=':', c='k')
        ax.text(
                0.95, 0.85, "$F_\\nu \propto \\nu^{%s \pm %s}$" %(m,em), 
                transform=ax.transAxes, horizontalalignment='right', fontsize=14)

        # the lines join at (144, 92)
        #ax.text(
        #        0, 0.9, "$(\\nu_p, F_p) = (144\,\mathrm{GHz}, 92\,\mathrm{mJy}$)", 
        #        transform=ax.transAxes, horizontalalignment='left', fontsize=14)

    if quad:
        # quadratic function to find peak
        nuc = np.logical_and(nu>87.5, nu<107.5)
        out = np.polyfit(nu[nuc], flux[nuc], deg=2)
        xlab = np.linspace(87.5, 107.5)
        ylab = out[0]*xlab**2 + out[1]*xlab + out[2]
        ax.plot(xlab, ylab, c='k', ls='--')

    if day == 14:
        day = '13 \& 14'
    ax.text(
            0.9, 0.1, "Day %s" %day, 
            transform=ax.transAxes, horizontalalignment='right', 
            fontsize=14)

    ax.set_ylim(0.1, 200)
    ax.tick_params(axis='both', labelsize=14)
    ax.tick_params(axis='both', labelsize=14)


if __name__=="__main__":
    fig, axarr = plt.subplots(
            3, 1, figsize=(8,8), sharex=True, sharey=True)

    # top panel: spectrum on Day 10
    ax = axarr[0]
    tel, freq, day, flux, eflux_form, eflux_sys = get_data_all()
    choose = day == 10
    islim = np.array([False]*sum(choose))
    islim[freq[choose]==5.5] = True
    plot_day(ax,10,freq[choose],flux[choose],islim)

    # middle panel: spectrum on Day 14
    ax = axarr[1]
    tel, freq, day, flux, eflux_form, eflux_sys = get_data_all()
    choose_cm = np.logical_and(day==13, tel=='ATCA')
    choose_mm = np.logical_and(
            day==14, 
            np.logical_or(tel=='SMA', tel=='ALMA'))
    choose = np.logical_or(choose_cm, choose_mm)
    islim = np.array([False]*sum(choose))
    plot_day(ax,14,freq[choose],flux[choose],islim)

    # bottom panel: spectrum on Day 22
    ax = axarr[2]
    freq,flux = get_spectrum(22)
    islim = np.array([False]*(len(freq)))
    plot_day(ax, 22, freq, flux, islim)
    ax.set_xlabel("Frequency [GHz]", fontsize=16)
    ax.text(
            0.1, 0.93, "$\\nu_p \\approx 100, F_p \\approx 94$", 
            fontsize=14, transform=ax.transAxes, 
            horizontalalignment='left', verticalalignment='top') 

    # Zoomed-in window
    axins = inset_axes(
            ax, 2, 1, loc=1,
            bbox_to_anchor=(0.4,0.8),
            bbox_transform=ax.transAxes)
    plot_day(
            axins, 22, freq, flux, islim,
            formatting=False, fit_peak=False, quad=True)
    axins.set_xlim(88, 106)
    axins.set_ylim(90, 95)
    axins.tick_params(axis='both', labelsize=12)
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    axins.xaxis.set_major_locator(MaxNLocator(nbins=2, prune='lower'))
    axins.yaxis.set_major_locator(MaxNLocator(nbins=2, prune='lower'))



    plt.show()
    #plt.savefig("spec.png")
