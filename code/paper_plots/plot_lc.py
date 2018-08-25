import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib
import matplotlib.gridspec as gridspec
import numpy as np
from astropy.table import Table
from xray_lc import get_xrt


def sma(ax):
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
    eflux_sys = np.array([0.1*f for f in flux])
    eflux_form = np.array(
            [float(val.split("pm")[1][0:-1]) for val in flux_raw])
    eflux = np.sqrt(eflux_sys**2 + eflux_form**2)

    #choose = np.logical_and(freq > 230, freq < 245)
    choose = freq == 231.5

    plot_flux = flux[choose]
    plot_eflux = eflux[choose]

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

    ax.set_ylim(4,70)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel("$F_{\\nu}$ [mJy]", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.legend(fontsize=12, loc='lower left')
    ax.set_yticks([10,50])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())



def atca(ax):
    # Get the X-ray data
    x,y,yerr = get_xrt() # in erg/cm^2/s
    f = (y / 10**(17)) * 1e23 * 1e3
    fact = 100
    ax.scatter(
            x, fact*f, c='grey', s=10)
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
    eflux_sys = np.array(
            [0.1*f for f in flux])
    eflux_form = np.array(
            [float(val.split("pm")[1][0:-1]) for val in flux_raw])
    eflux = np.sqrt(eflux_sys**2 + eflux_form**2)
    eflux[eflux_form < 0] = -99 # non-detections

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
    ax.locator_params(axis='y', nbins=2)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.set_yscale('log')
    ax.legend(fontsize=12, loc='lower center', ncol=3)
    ax.set_yticks([0.1,1])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())



if __name__=="__main__":
    """ Plot ATCA and SMA light curves """
    fig = plt.subplots(figsize=(6,6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3,2], hspace=0.0)
    gs.update(left=0.15, right=0.9)

    atca_ax = plt.subplot(gs[1])
    sma_ax = plt.subplot(gs[0], sharex=atca_ax)
    atca(atca_ax)
    sma(sma_ax)
    plt.setp(sma_ax.get_xticklabels(), visible=False)
    #plt.tight_layout()

    #plt.savefig("lc.png")
    plt.show()
