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



def plot_day(ax,day,nu,flux,tel,islim,formatting=True,fit_peak=True,quad=False,showlegend=False):
    """ Run this for one day """
    choose_tel = tel == 'SMA'
    choose = np.logical_and(choose_tel, ~islim)
    ax.errorbar(
            nu[choose], flux[choose], yerr=0.15*flux[choose],
            fmt='o', label="SMA", c='k', zorder=10)
    choose = np.logical_and(choose_tel, islim)
    ax.scatter(
            nu[choose], flux[choose], marker='v', c='k')

    choose_tel = tel == 'ALMA'
    choose = np.logical_and(choose_tel, ~islim)
    ax.scatter(
            nu[choose], flux[choose], 
            marker='o', label="ALMA", facecolor='#f98e09',
            edgecolor='k', zorder=5)
    choose = np.logical_and(choose_tel, islim)
    ax.scatter(
            nu[choose], flux[choose], marker='v', c='#f98e09')

    choose_tel = tel == 'ATCA'
    choose = np.logical_and(choose_tel, ~islim)
    ax.scatter(
            nu[choose], flux[choose], 
            marker='s', label="ATCA", facecolor='#57106e',
            edgecolor='#57106e', zorder=5)
    choose = np.logical_and(choose_tel, islim)
    ax.scatter(
            nu[choose], flux[choose], marker='v', c='#57106e')

    if formatting:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylabel("Sp. Flux Dens. [mJy]", fontsize=14)
        ax.text(
                0.35, 0.45, "$F_\\nu \propto \\nu^{2.5}$", 
                transform=ax.transAxes, horizontalalignment='right', 
                fontsize=14)
    if fit_peak:
        choose = np.logical_and(np.logical_and(nu>12, nu<70), ~islim)
        b,berr = fit_self_abs(nu[choose],flux[choose],0.1*flux[choose])
        b = b[0]
        xfit = np.linspace(5, 200)
        yfit = 10**(2.5*np.log10(xfit)+b)
        ax.plot(xfit, yfit, ls=':', c='k')

        # fit and plot a nu^something line
        xfit = np.linspace(50, 1000)
        if day == 10:
            m = -1.86
            em = 0.03
            b = 5.9824
        elif day == 14:
            m = -0.75
            em = 0.01
            b = 3.36490402919
        elif day == 22:
            m = -1.06
            em = 0.01
            b = 4.20478080158
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

    if day == 10:
        day = '=10.2--10.5'
    if day == 14:
        day = '=13.36--14.38'
    if day == 22:
        day = '$\\approx 22$'
    ax.text(
            0.9, 0.1, "$\Delta t\,$[days]%s" %day, 
            transform=ax.transAxes, horizontalalignment='right', 
            fontsize=14)

    if showlegend:
        ax.legend(loc='upper left', fontsize=12)

    ax.set_ylim(0.1, 200)
    ax.tick_params(axis='both', labelsize=14)
    ax.tick_params(axis='both', labelsize=14)


if __name__=="__main__":
    fig, axarr = plt.subplots(
            3, 1, figsize=(8,8), sharex=True, sharey=True, dpi=100)

    # spectrum on Day 10
    # 10.26 for SMA, 10.48 for ATCA
    ax = axarr[0]
    tel, freq, day, flux, eflux_form, eflux_sys = get_data_all()
    choose = np.logical_and(day > 10.2, day < 10.5)
    islim = np.array([False]*sum(choose))
    islim[freq[choose]==5.5] = True
    plot_day(ax,10,freq[choose],flux[choose],tel[choose],islim, showlegend=True)

    # spectrum on Day 13/14
    # 14.36 for SMA, 13.47 for ATCA, 14.14 for ALMA
    ax = axarr[1]
    tel, freq, day, flux, eflux_form, eflux_sys = get_data_all()
    choose = np.logical_and(day < 14.38, day > 13.36)
    islim = np.array([False]*sum(choose))
    plot_day(ax,14,freq[choose],flux[choose],tel[choose],islim)

    # bottom panel: spectrum on Day 22
    # 22.02 to 22.04 for ALMA
    ax = axarr[2]
    choose = np.logical_and(day < 22.05, day > 22.01) 
    nu = freq[choose]
    f = flux[choose]
    t = tel[choose]
    # interpolate the SMA data between Days 20 and 24
    choose_sma = np.logical_and.reduce((
        day>=20.28, day<=24.39, tel=='SMA', freq==215.5))
    nu = np.append(nu, 215.5)
    f = np.append(f, np.interp(22, day[choose_sma], flux[choose_sma]))
    t = np.append(t, 'SMA')
    choose_sma = np.logical_and.reduce((
        day>=20.28, day<=24.39, tel=='SMA', freq==231.5))
    nu = np.append(nu, 231.5)
    f = np.append(f, np.interp(22, day[choose_sma], flux[choose_sma]))
    t = np.append(t, 'SMA')
    # interpolate the ATCA 34 GHz measurement
    choose_atca = np.logical_and.reduce((
        tel=='ATCA', freq==34))
    nu = np.append(nu, 34)
    f = np.append(f, np.interp(22, day[choose_atca], flux[choose_atca]))
    t = np.append(t, 'ATCA')

    islim = np.array([False]*(len(nu)))
    plot_day(ax, 22, nu, f, t, islim)
    ax.set_xlabel("Frequency [GHz]", fontsize=14)
    ax.scatter(
            671, 31.5, marker='*', s=100, 
            facecolor='white', edgecolor='black')
    ax.text(
            750, 25, 'Day 24', fontsize=10)
    ax.text(
            0.1, 0.93, "$\\nu_p \\approx 100, F_p \\approx 94$", 
            fontsize=14, transform=ax.transAxes, 
            horizontalalignment='left', verticalalignment='top') 

    # Zoomed-in window
    # axins = inset_axes(
    #         ax, 2, 1, loc=1,
    #         bbox_to_anchor=(0.4,0.8),
    #         bbox_transform=ax.transAxes)
    # plot_day(
    #         axins, 22, nu, f, t, islim,
    #         formatting=False, fit_peak=False, quad=True)
    # axins.set_xlim(88, 106)
    # axins.set_ylim(90, 95)
    # axins.tick_params(axis='both', labelsize=12)
    # mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    # axins.xaxis.set_major_locator(MaxNLocator(nbins=2, prune='lower'))
    # axins.yaxis.set_major_locator(MaxNLocator(nbins=2, prune='lower'))

    plt.show()
    #plt.savefig("spec.eps", format='eps', bbox_inches='tight', pad_inches=0.1)
