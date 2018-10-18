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

    a, b, c = sma_lc()
    dt, f, ef = b # 230 GHz light curve
    ef_comb = np.sqrt(ef**2 + (0.15*f)**2)
    ax.scatter(
            dt[:-1], f[:-1],
            marker='s', c='k',
            label="230.6--234.6 GHz")
    ax.errorbar(
            dt[:-1], f[:-1], ef_comb[:-1], 
            fmt='s', c='#000004', lw=1.5) # black in inferno
    ax.plot(
            dt[:-1], f[:-1], linestyle='-', c='k', lw=1.5)
    plt.axvline(x=16)

    dt, f, ef = c # 345 GHz light curve
    ef_comb = np.sqrt(ef**2 + (0.15*f)**2)
    ax.errorbar(
            dt, f, ef_comb, 
            fmt='*', c='#f98e09', lw=0.5, alpha=0.8, ms=11,
            label="341.5--349 GHz")
    ax.plot(
            dt, f, linestyle='-', c='#f98e09', lw=1.5)

    # non-detection
    # what is the 3-sigma of the non-detection?
    # x + 3-sigma
    # for 231.5 GHz, it's X=0.62, 3-sig = 3*0.64
    val = 0.62
    sig = 0.64
    ax.scatter(76.27, val+3*sig, marker="_", c='k', s=200)
    ax.scatter(76.27, val, marker="_", c='k', s=100)
    ax.scatter(76.27, val, marker=".", c='k', s=30)
    ax.arrow(
            76.27, val+3*sig, 0, -3*sig, length_includes_head=True, 
            head_width=5, head_length=0.1, fc='k')
    

    # for 351.0, the non-detection is -0.32 +/- 1.76
    val = -0.32
    sig = 1.76
    ax.scatter(76.27, val+3*sig, marker="_", c='#f98e09', s=200)
    ax.scatter(76.27, val, marker="_", c='#f98e09', s=100)
    ax.scatter(76.27, val, marker=".", c='#f98e09', s=30)
    ax.arrow(
            76.27, val+3*sig, 0, -3*sig, length_includes_head=True, 
            head_width=5, head_length=0.1, fc='#f98e09', ec='#f98e09')

    # Do 34 GHz
    choose = freq == 34
    # dont include the last point
    xfit = days[choose][0:-1]
    yfit = flux[choose][0:-1]
    eyfit = eflux[choose][0:-1]
    # fit a nu^2 power law
    out = np.polyfit(np.log10(xfit), np.log10(yfit), deg=1)
    x = np.linspace(min(xfit), max(xfit))
    y = 10**(out[0]*np.log10(x)+out[1])
    ax.errorbar(
            days[choose], flux[choose], yerr=eflux[choose],
            fmt='o', c='#57106e', lw=1.5, ms=7,
            label="34 GHz")
    ax.plot(x, y, ls='--', c='#57106e')
    ax.text(17, 9.5, r'$f_\nu \propto t^2$', fontsize=12)

    # for SMA proposal
    #ax.axvline(x=5, color='#8da0cb', alpha=0.5, lw=7)
    #ax.axvline(x=12, color='#8da0cb', alpha=0.5, lw=7)

    # cross hatches for the Day 10, Day 13/14, and Day 22 measurements
    ax.text(10, 70, 'S', fontsize=12)
    ax.text(14, 70, 'S', fontsize=12)
    ax.text(22, 70, 'S', fontsize=12)

    ax.set_ylim(0.5,100)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.legend(fontsize=12, loc='lower left',ncol=3)
    ax.set_yticks([0.5, 10, 50])
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
    print(dt)
    ax.errorbar(
            dt, f1*1E12, yerr=ef1*1E12, fmt='o', label="3-10 keV",
            mfc='#bc3754', mec='k', c='k', ms=msize)
    ax.errorbar(
            dt, f2*1E12, yerr=ef2*1E12,
            fmt='v', mfc='white', mec='black',
            label="10-20 keV", c='k', ms=msize)
    ax.errorbar(
            dt, f3*1E12, yerr=ef3*1E12,
            fmt='^', label="20-40 keV",
             mfc='#f98e09', mec='black', c='k', ms=msize)
    ax.errorbar(
            dt, f4*1E12, yerr=ef4*1E12, fmt='s',
            mfc='#57106e', mec='black',
            label="40-80 keV", c='k', ms=msize)

    ax.set_xlabel("Days since MJD 58285.0 (2018 June 16 UT)", fontsize=16)
    ax.set_xlim(2.5,100)
    ax.set_ylabel("$F_{X}$ [$10^{-12}$ erg/cm${}^2$/s]", fontsize=16)
    ax.locator_params(axis='y', nbins=2)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)
    ax.set_yscale('log')
    ax.legend(fontsize=11, loc='lower left', ncol=3)
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
