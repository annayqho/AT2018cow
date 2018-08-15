import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.gridspec as gridspec
import numpy as np
from astropy.table import Table
from xray_lc import get_xray


def spindex(ax, d, nulow, nuhigh, flow, eflow, fhigh, efhigh):
    alpha = np.log10( fhigh / flow ) / np.log10( nulow / nuhigh )

    # estimate the uncertainty
    a = (1 / np.log10(nulow / nuhigh)) 
    b = flow/fhigh
    c = np.log10(np.e)
    dadf_hi = a*b*c* (1/flow) * efhigh
    dadf_low = a*b*c* (-fhigh/flow**2) * eflow
    ealpha = np.sqrt(dadf_hi**2 + dadf_low**2)

    #ax.scatter(d, alpha, c='k')
    ax.errorbar(d, alpha, yerr=ealpha, c='k', fmt='.')
    ax.plot(d, alpha, c='k', lw=1.0)
    ax.set_ylabel("$\\alpha$", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)

    return alpha, ealpha


def sma(ax, legend, plot230=True, plot340=True, unit='mjy'):
    data_dir = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
    dat = Table.read(
        "%s/radio_lc.dat" %data_dir, delimiter="&", format='ascii.no_header')
    tel = np.array(dat['col2'])
    choose = np.logical_or(tel == 'SMA', tel == 'ATCA')

    days = np.array(dat['col1'][choose])
    freq = np.array(dat['col3'][choose]).astype(float)
    flux_raw = np.array(dat['col4'][choose])
    flux = np.array(
            [float(val.split("pm")[0][1:]) for val in flux_raw])
    eflux = np.array(
            [float(val.split("pm")[1][0:-1]) for val in flux_raw])

    #choose = np.logical_and(freq > 230, freq < 245)
    choose = freq == 231.5

    if plot230:
        if unit=='cgs':
            plot_flux = flux[choose] * 1e-3 * 1e-23 # cgs units from mJy
            plot_eflux = eflux[choose] * 1e-3 * 1e-23
        elif unit=='mjy':
            plot_flux = flux[choose]
            plot_eflux = eflux[choose]
        elif unit=='nufnu':
            plot_flux = flux[choose] * 1e-3 * 1e-23 * freq[choose] * 1e9
            plot_eflux = eflux[choose] * 1e-3 * 1e-23 * freq[choose] * 1e9
        else:
            print("error: I don't recognize this unit")
        # just for the label
        ax.scatter(
                days[choose], plot_flux, 
                marker='s', c='k',
                label="231.5 GHz")
        ax.errorbar(
                days[choose], plot_flux, plot_eflux, 
                fmt='s', c='k', lw=1.5)
        ax.plot(
                days[choose], plot_flux, linestyle='-', c='k', lw=1.5)


    choose = np.logical_and(freq >= 341.5, freq <= 349)

    if plot340:
        #ax.scatter(
        #        days[choose], flux[choose], 
        #        marker='*', c='grey', alpha=0.8, s=13,
        ax.errorbar(
                days[choose], flux[choose], eflux[choose],
                fmt='*', c='grey', lw=0.5, alpha=0.8, ms=13,
                label="341.5--349 GHz")
        ax.plot(
                days[choose], flux[choose], linestyle='-', c='grey', lw=1.5)

    # Do 34 GHz
    choose = freq == 34
    xfit = days[choose]
    yfit = flux[choose]
    eyfit = eflux[choose]
    # fit a nu^2 power law
    out = np.polyfit(np.log10(xfit), np.log10(yfit), deg=1)
    x = np.linspace(min(xfit), max(xfit))
    y = 10**(out[0]*np.log10(x)+out[1])
    ax.errorbar(
            xfit, yfit, yerr=eyfit,
            fmt='o', c='k', lw=1.5, ms=7,
            label="34 GHz")
    ax.plot(x, y, ls='--', c='k')

    if np.logical_or(plot230, plot340):
        ax.set_ylim(4,70)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_ylabel("$F_{\\nu}$ [mJy]", fontsize=16)
        ax.yaxis.set_tick_params(labelsize=14)

        #ax.text(0.05, 0.9, "SMA", transform=ax.transAxes, 
        #        fontsize=, verticalalignment='top')

        if legend:
            ax.legend(fontsize=12, loc='lower left')

def atca(ax):
    # Get the X-ray data
    x,y,yerr = get_xray() # in erg/cm^2/s
    f = (y / 10**(17)) * 1e23 * 1e3
    fact = 100
    ax.plot(
            x, fact*f, lw=0.5, c='grey', 
            label="0.3-10 keV($\\times %s$)" %fact)

    ax.set_xlabel("Time Since June 16 UT [d]", fontsize=16)
    data_dir = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
    dat = Table.read(
        "%s/radio_lc.dat" %data_dir, delimiter="&", format='ascii.no_header')
    tel = np.array(dat['col2'])
    choose = tel == 'ATCA'

    days = np.array(dat['col1'][choose])
    freq = np.array(dat['col3'][choose]).astype(float)
    flux_raw = np.array(dat['col4'][choose])
    flux = np.array(
            [float(val.split("pm")[0][1:]) for val in flux_raw])
    eflux = np.array(
            [float(val.split("pm")[1][0:-1]) for val in flux_raw])

    eflux[eflux < 0] = -99 # non-detections

    choose_freq = freq==5.5

    # for detections
    det = eflux > 0
    choose = np.logical_and(choose_freq, det)

    ax.errorbar(
            days[choose], flux[choose], eflux[choose], fmt='s', c='grey',
            label = "5.5 GHz")
    xfit = days[choose]
    yfit = flux[choose]
    eyfit = eflux[choose]
    # fit a nu^2 power law
    out = np.polyfit(np.log10(xfit), np.log10(yfit), deg=1)
    x = np.linspace(min(xfit), max(xfit))
    y = 10**(out[0]*np.log10(x)+out[1])
    ax.plot(x, y, ls='--', c='grey')

    ax.set_xlabel("Time Since June 16 UT [d]", fontsize=16)

    # for non-detections
    choose = np.logical_and(choose_freq, ~det)

    ax.scatter(
            days[choose], flux[choose], c='grey', marker='v')

    choose = freq == 9

    xfit = days[choose]
    yfit = flux[choose]
    eyfit = eflux[choose]
    # fit a nu^2 power law
    out = np.polyfit(np.log10(xfit), np.log10(yfit), deg=1)
    x = np.linspace(min(xfit), max(xfit))
    y = 10**(out[0]*np.log10(x)+out[1])
    ax.errorbar(
            xfit, yfit, yerr=eyfit,
            fmt='o', c='k', lw=1.5, ms=7,
            label="9 GHz")
    ax.plot(x, y, ls='--', c='k')

    ax.set_xlim(3,50)
    ax.set_ylabel("$F_{\\nu}$ [mJy]", fontsize=16)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)

    ax.set_yscale('log')
    #ax.text(
    #        0.05, 0.9, "ATCA", transform=ax.transAxes, fontsize=14,
    #        verticalalignment='top')

    ax.legend(fontsize=12, loc='lower center', ncol=3)


if __name__=="__main__":
    """ Plot ATCA and SMA light curves """
    fig = plt.figure(figsize=(6,6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3,2], hspace=0.0)
    gs.update(left=0.15, right=0.9)

    atca_ax = plt.subplot(gs[1])
    sma_ax = plt.subplot(gs[0], sharex=atca_ax)
    atca(atca_ax)
    sma(sma_ax, legend=True)
    plt.setp(sma_ax.get_xticklabels(), visible=False)

    plt.savefig("lc.png")
    #plt.show()
