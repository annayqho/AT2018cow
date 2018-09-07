import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib
import matplotlib.gridspec as gridspec
import numpy as np
from astropy.table import Table
from astropy.time import Time
import sys
sys.path.append('/Users/annaho/Dropbox/Projects/Research/AT2018cow/code')
from xray_lc import get_xrt, get_nustar
from scale_fluxes import sma_lc


def sma(ax):
    data_dir = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
    dat = Table.read(
        "%s/radio_lc.dat" %data_dir, delimiter="&", 
        format='ascii.no_header')
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

    a, b = sma_lc()
    dt, f, ef = b
    ef_comb = np.sqrt(ef**2 + (0.15*f)**2)
    ax.scatter(
            dt, f,
            marker='s', c='k',
            label="230.6--234.6 GHz")
    ax.errorbar(
            dt, f, ef_comb, 
            fmt='s', c='k', lw=1.5)
    ax.plot(
            dt, f, linestyle='-', c='k', lw=1.5)

    choose = np.logical_and(freq >= 341.5, freq <= 349)

    ax.errorbar(
            days[choose], flux[choose], eflux[choose],
            fmt='*', c='#66c2a5', lw=0.5, alpha=0.8, ms=13,
            label="341.5--349 GHz")
    ax.plot(
            days[choose], flux[choose], linestyle='-', c='#66c2a5', lw=1.5)

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
            fmt='o', c='#8da0cb', lw=1.5, ms=7,
            label="34 GHz")
    ax.plot(x, y, ls='--', c='#8da0cb')

    # for SMA proposal
    #ax.axvline(x=5, color='#8da0cb', alpha=0.5, lw=7)
    #ax.axvline(x=12, color='#8da0cb', alpha=0.5, lw=7)

    ax.set_ylim(4,70)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.legend(fontsize=12, loc='lower left')
    ax.set_yticks([10,50])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())



def xray(ax):
    # Get the X-ray data
    x,y,yerr = get_xrt() # in erg/cm^2/s
    ax.errorbar(
            x, y*1E12, yerr=yerr*1E12, c='k', ms=3, label=None, fmt='.')
    ax.plot(
            x, y*1E12, lw=0.5, c='k', 
            label="XRT 0.3-10 keV")

    msize=7
    dt, f1, ef1, f2, ef2, f3, ef3, f4, ef4 = get_nustar()
    ax.errorbar(
            dt, f1*1E12, yerr=ef1*1E12, fmt='o', label="3-10 keV",
            mfc='white', mec='k', c='k', ms=msize)
    ax.errorbar(
            dt, f2*1E12, yerr=ef2*1E12,
            fmt='v', mfc='#1b9e77', mec='black',
            label="10-20 keV", c='k', ms=msize)
    ax.errorbar(
            dt, f3*1E12, yerr=ef3*1E12,
            fmt='^', label="20-40 keV",
             mfc='#d95f02', mec='black', c='k', ms=msize)
    ax.errorbar(
            dt, f4*1E12, yerr=ef4*1E12, fmt='s',
            mfc='#7570b3', mec='black',
            label="40-80 keV", c='k', ms=msize)

    ax.set_xlabel("Time Since June 16 UT (MJD 58285) [d]", fontsize=16)
    ax.set_xlim(2.5,80)
    ax.set_ylabel("$F_{X}$ [$10^{-12}$ erg/cm${}^2$/s]", fontsize=16)
    ax.locator_params(axis='y', nbins=2)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.set_yscale('log')
    ax.legend(fontsize=11, loc='lower left', ncol=2)
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())



if __name__=="__main__":
    """ Plot radio and X-ray light curves """
    fig = plt.subplots(figsize=(8,6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[5,4], hspace=0.0)
    gs.update(left=0.15, right=0.9)

    xray_ax = plt.subplot(gs[1])
    sma_ax = plt.subplot(gs[0], sharex=xray_ax)
    xray(xray_ax)
    sma(sma_ax)
    plt.setp(sma_ax.get_xticklabels(), visible=False)
    #plt.tight_layout()

    #plt.savefig("lc.png")
    plt.show()
