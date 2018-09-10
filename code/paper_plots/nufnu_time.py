""" 
Plot of luminosity over time
"""


import matplotlib
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/code")
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/data/radio_compilations/Zauderer2011")
from astropy.cosmology import Planck15
from get_radio import *
from scale_fluxes import sma_lc
from read_table import *


def plot_limits(ax, x, y, ratiox, ratioy, col):
    """ Plot two arrows from the point """
    ax.annotate('', xy=(x*ratiox, y), xytext=(x, y),
            arrowprops=dict(
                facecolor=col, headwidth=10, width=1, headlength=7))
    ax.annotate('', xy=(x, y*ratioy), xytext=(x, y),
            arrowprops=dict(
                facecolor=col, headwidth=10, width=1, headlength=7))


def plot_line(ax, d, nu, t, f, name, label, col, legend=False):
    """ Plot a line
    If nu > 90 GHz, make it on the left axis
    If nu < 10 GHz, make it on the right axis
    
    Parameters
    ----------
    f: flux in mJy
    name: name of the source
    label: label to use as the legend (which also determines the col)
    """
    lum = nu * f * 1e-23 * 1e-3 * 4 * np.pi * d**2
    fs = 11
    nsize = 10 # normal size for points
    if name=='AT2018cow':
        marker='*'
        fcol = col
        s=70
    else:
        if label=='SN':
            marker='o'
            s=nsize
            fcol = col # fill color
        elif label=='GRB':
            marker='o'
            fcol = 'white' # unfilled
            s=nsize
        elif label=='Rel. SN':
            marker='s'
            fcol = col 
            s=nsize
        elif label=='Rel. TDE':
            marker='s'
            fcol = 'white' #unfilled
            s=nsize
    ax.scatter(
            t, lum, facecolor=fcol, edgecolor=col, 
            marker=marker, s=s)
    if legend:
        ax.plot(t, lum, c=col, ls='-', label=label)
    else:
        ax.plot(t, lum, c=col, ls='-', label=None)
    return lum


def plot_points(ax, d, nu, t, f, marker, name=None):
    """ Plot set of two points """
    lums = []
    for ii,nuval in enumerate(nu):
        if nuval > 90E9:
            lum = plot_point(ax, d, nuval, t[ii], f[ii], marker, name=name)
            lums.append(lum)
        else:
            lum = plot_point(ax, d, nuval, t[ii], f[ii], marker)
            lums.append(lum)
    ax.plot(
        t, lums, ls='--', c='k')
    return lums


def at2018cow(ax, col, legend):
    """ 231.5 GHz light curve and 9 GHz light curve """
    d = Planck15.luminosity_distance(z=0.014).cgs.value

    # high frequency
    a, b = sma_lc()
    dt, f, ef = b
    ef_comb = np.sqrt(ef**2 + (0.15*f)**2)
    nu = 231.5E9
    lum = plot_line(
            ax[0], d, nu, dt, f, 
            'AT2018cow', None, col, legend)
    ax[0].text(dt[0]/1.05, lum[-1], 'AT2018cow', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')

    # low frequency
    nu = 34E9
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
    choose = freq == 34
    lum = plot_line(
            ax[1], d, nu, days[choose], flux[choose], 
            'AT2018cow', None, col, legend)
    ax[1].text(days[choose][0]/1.1, lum[-1], 'AT2018cow', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


def maxi(ax):
    """ data from Dillon """
    d = 4.62E26
    t = 4*365
    nu = 6E9
    f = 6.16
    lum = plot_point(ax, d, nu, t, f, '*')
    plt.text(t/1.2, lum, 'MAXI 140814A', fontsize=11,
            horizontalalignment='right', verticalalignment='center')


def tde(ax, col, legend):
    """  Plot the 225 GHz light curve from the SMA
    
    Plot the 4.9 GHz light curve from the VLA
    """
    z = 0.354
    d = Planck15.luminosity_distance(z=z).cgs.value
    nu, dt, f, ef, islim = zauderer()
    t = dt/(1+z)

    nu_plt = 225E9
    choose = np.logical_and(~islim, nu == nu_plt/1E9)
    plot_line(
            ax[0], d, nu_plt, t[choose], f[choose], 
            'SwiftJ1644+57', 'Rel. TDE', col, legend)

    nu_plt = 4.9E9
    choose = np.logical_and(~islim, nu == nu_plt/1E9)
    plot_line(
            ax[1], d, nu_plt, t[choose], f[choose], 
            'SwiftJ1644+57', 'Rel. TDE', col, legend)


def sn2003L(ax, col, legend):
    """ Soderberg et al
    Values at 8.5 GHz """
    d = 2.8432575937224894e+26
    nu_plt = 8.5E9
    nu, dt, f, ef = read_2003L()
    choose = nu == nu_plt
    plot_line(
            ax[1], d, nu_plt, dt[choose], 1E-3*f[choose], 
            'SN2003L', 'SN', col, legend)
    

def sn1979c(ax, col, legend):
    """ Weiler 1986 and Weiler 1991
    This is a IIL 
    
    Too old to have the raw table on arXiv
    """
    d = 5.341805643483106e+25
    nu = 1.4E9 # 20cm
    t = np.array(
            [437,594,631,663,679,684,727,747,786,822,839,876,882,
             914,937,973,995,
             1026,1071,1091,1127,1156,1168,1212,1243,1277,1314,1358,1390,
             1415,1435,1466,1513,1565,1600,1634,1659,1698,1714,1750,1771,
             1931,2027])
    flux = np.array(
            [0.2,2.1,2.5,2.7,2.8,2.8,4.4,4.8,6.0,7.1,7.1,7.6,8.6,
             9.8,6.5,8.6,9.5,
             10.2,10.8,10.3,10.4,12.2,10.1,10.2,11.5,11.2,13.0,11.3,10.2,
             9.6,11.2,13.2,11.1,9.1,8.5,9.1,8.8,10.1,9.7,9.1,8.9,
             7.0,7.7])
    plot_line(ax[1], d, nu, t, flux, 'SN1979c', 'SN', col, legend)
    

def sn1993J(ax, col, legend):
    """ SN 1993J from Weiler et al. 
    This is the peak of the 99.4 GHz light curve
    There is also a 110 GHz point, but only one,
    so I can't get the peak.
    1.4 GHz / 20 cm, 4.9 GHz / 6 cm
    values come from the best-fit model,
    but by eye they are clearly pretty close
    """
    d = 1.1e25
    freq = 5E9
    nu, dt, f, ef, islim = read_1993J_low_freq()
    choose = np.logical_and(~islim, nu==freq)
    lum = plot_line(
            ax[1], d, freq, dt[choose], f[choose], 
            'SN1993J', 'SN', col, legend)
    ax[1].text(dt[choose][0]/1.05, lum[0], 'SN1993J', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')

    freq = 99.4E9
    nu, dt, f, ef, islim = read_1993J_high_freq()
    choose = np.logical_and(~islim, nu==freq)
    lum = plot_line(
            ax[0], d, freq, dt[choose], f[choose], 
            'SN1993J', 'SN', col, legend)
    ax[0].text(dt[choose][-1]*1.1, lum[-1], 'SN1993J', fontsize=11,
            verticalalignment='center',
            horizontalalignment='left')


def sn2011dh(ax, col, legend):
    """ SN 2011dh
    Horesh et al. 
    M51: d = 8.03 Mpc; expl date May 31.58
    The high-freq points are limits. The peak was earlier and more luminous.
    The low-freq points are limits. The peak was later and more luminous.
    """
    d = 2.5E25
    freq = 107E9
    dt, nu, f, ef, islim = read_2011dh()
    choose = np.logical_and(~islim, nu==freq)
    lum = plot_line(
            ax[0], d, nu[choose], dt[choose], f[choose], 
            'SN2011dh', 'SN', col, legend)
    ax[0].text(dt[choose][0]/1.05, lum[0], 'SN2011dh', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')

    freq = 8.5E9
    choose = np.logical_and(~islim, nu==freq)
    lum = plot_line(
            ax[1], d, nu[choose], dt[choose], f[choose], 
            'SN2011dh', 'SN', col, legend)
    ax[1].text(dt[choose][0]/1.05, lum[0], 'SN2011dh', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


def grb030329(ax, col, legend):
    """ Sheth et al. 2003 
    
    Explosion day was obviously 03/29
    """
    freq = 250E9
    z = 0.1686
    d = Planck15.luminosity_distance(z=z).cgs.value

    t = np.array([1, 4, 6, 7, 8, 16, 22]) / (1+z)
    f = np.array([49.2, 45.7, 41.6, 32, 25.5, 9.3, 5.2])
    plot_line(ax[0], d, freq, t, f, 'GRB030329', 'GRB', col, legend)

    freq = 8.5E9
    t = np.array(
            [0.58, 1.05, 2.65, 3.57, 4.76, 6.89, 7.68, 9.49, 11.90, 
                12.69, 14.87, 16.66, 18.72, 20.58, 25.70, 28.44, 31.51, 
                33.58, 36.52, 42.55, 44.55, 59.55, 66.53]) / (1+z)
    f = np.array(
            [3.50, 1.98, 8.50, 6.11, 9.68, 15.56, 12.55, 13.58, 17.70, 
                17.28, 19.15, 17.77, 15.92, 16.08, 15.34, 12.67, 13.55, 
                13.10, 10.64, 8.04, 8.68, 4.48, 4.92])
    plot_line(ax[1], d, freq, t, f, 'GRB030329', 'GRB', col, legend)


def grb130427A(ax, col, legend):
    """ Perley et al
    They have data from CARMA/PdBI at 90 GHz (3mm)
    But by the time they caught it, it was fading
    """
    z = 0.340
    d = Planck15.luminosity_distance(z=z).cgs.value

    freq = 93E9
    t = np.array([0.77, 1, 1.91, 2.8]) / (1+z)
    flux = np.array([3416, 2470, 1189, 807]) * 1E-3
    plot_line(ax[0], d, freq, t, flux, 'GRB130427A', 'GRB', col, legend)

    freq = 5.10E9
    t = np.array([0.677, 2.04, 4.75, 9.71, 17.95, 63.78, 128.34]) / (1+z)
    f = np.array([1290, 1760, 648, 454, 263, 151, 86]) * 1E-3
    plot_line(ax[1], d, freq, t, f, 'GRB130427A', 'GRB', col, legend)


def sn2007bg(ax, col, legend):
    """ Salas et al. 2013
    Peak is resolved for 4.86, 8.46 GHz
    Let's choose the one where it's more clearly resolved
    It's a bit weird because there are two peaks, but let's
    choose the first one because it neglects subsequent interaction
    """
    nu = 8.46E9
    d = Planck15.luminosity_distance(z=0.0346).cgs.value
    t = np.array(
            [13.8, 19.2, 26.1, 30.9, 41.3, 55.9, 66.8, 81.8, 98.8, 124, 
                144, 159.8, 189.9, 214.9, 250.9, 286.8, 314.8, 368.8, 
                386.8, 419.9, 566.9, 623.8, 720.8, 775.8, 863.8])
    f = np.array(
            [480, 753, 804, 728, 1257, 1490, 1390, 1325, 1131, 957, 
                621, 316, 379, 404, 783, 1669, 2097, 2200, 
                2852, 3344, 3897, 3891, 3842, 3641, 3408]) * 1E-3
    plot_line(ax[1], d, nu, t, f, 'SN2007bg', 'SN', col, legend)


def sn2003bg(ax, col, legend):
    """ Soderberg et al.
    Peak is resolved for 22.5, 15, 8.46, 4.86, 1.43
    Again, there are two peaks...
    Let's choose the first peak, 8.46
    """
    nu = 8.46E9
    d = 6.056450393620008e+25

    t = np.array(
                [10, 12, 23, 35, 48, 58, 63, 73, 85, 91, 115, 129,
                132, 142, 157, 161, 181, 201, 214, 227, 242, 255,
                266, 285, 300, 326, 337, 351, 368, 405, 410, 424,
                434, 435, 493, 533, 632, 702, 756, 820, 902, 978])
    f = np.array(
                [2.51, 3.86, 12.19, 24.72, 40.34, 51.72, 49.64, 46.20,
                38.638, 33.85, 45.74, 53.94, 54.27, 54.83, 48.43,
                47.43, 35.76, 31.35, 28.67, 27.38, 24.57, 22.30,
                21.67, 21.31, 20.88, 20.33, 19.85, 18.84, 17.14,
                14.61, 14.49, 14.16, 13.25, 13.08, 10.04, 8.92,
                6.23, 6.18, 4.62, 3.93, 4.69, 4.48])
    plot_line(ax[1], d, nu, t, f, 'SN2003bg', 'SN', col, legend)


def sn2009bb(ax, col, legend):
    """ expl date Mar 19 """
    nu = 8.46E9
    d = 1.237517263280789e+26
    t_apr = 11 + np.array([5.2, 8.2, 13.2, 15.1, 23.2, 29.1])
    t_may = 11 + 30 + np.array([3.1, 10.1, 13, 20.1, 27])
    t_jun = 11 + 30 + 31 + np.array([6, 17, 26, 18.9])
    t = np.hstack((t_apr, t_may, t_jun))
    flux = np.array([24.681, 17.568, 16.349, 13.812, 8.881,
        7.714, 8.482, 6.824, 6.327, 3.294, 4.204, 3.203, 2.392,
        1.903, 1.032])
    lum = plot_line(ax[1], d, nu, t, flux, 'SN2009bb', 'Rel. SN', col, legend)
    ax[1].text(t[0]/1.05, lum[0], '2009bb', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


def sn1998bw(ax, col, legend):
    """ SN 1998bw
    
    This is a bit complicated because there are two peaks in the light curve,
    but I am using a frequency that has a main peak rather than a frequecy
    with two clear distinct peaks...
    """
    d = 1.17E26 # cm
    nu = 150E9
    t = np.array([12.4])
    f = np.array([39])
    lum = plot_line(ax[0], d, nu, t, f, 'SN1998bw', 'Rel. SN', col, legend)
    nu = 2.3E9
    t = np.array([11.7, 14.6, 15.7, 16.5, 17.8, 19.7, 21.6, 23.6, 25.9, 26.8, 28.8, 30.0, 32.9, 34.7, 36.8, 38.8, 40.0, 45.7, 51.7, 57.7, 64.7, 67.7, 80.5])
    f = np.array([19.7, 22.3, 23.5, 23.9, 25.1, 25.3, 20.9, 22.9, 28.0, 28.7, 31.1, 31.3, 27.3, 33.5, 31.8, 31, 31.3, 26.8, 23.1, 18.5, 15.6, 15.6, 9.6])
    ax[0].text(t[0], lum[0]/1.3, '1998bw', fontsize=11,
            verticalalignment='top',
            horizontalalignment='center')
    lum = plot_line(ax[1], d, nu, t, f, 'SN1998bw', 'Rel. SN', col, legend)
    ax[1].text(t[0]/1.05, lum[0], '1998bw', fontsize=11,
            verticalalignment='center',
            horizontalalignment='right')


def othersn(ax):
    """ 
    SN 2013ak
    most luminous events are typically SNe IIn (Chevalier 2006)
    SN Ib/c: Berger 2003, Chevalier & Fransson 2006
    2008D
    """


if __name__=="__main__":
    fig, axarr = plt.subplots(1, 2, figsize=(10,5), sharex=True, sharey=True)
    props = dict(boxstyle='round', facecolor='white')

    # nice colors from matplotlib colormap "Paired"
    #a6cee3 light pastel blue
    #1f78b4 a more vibrant normal blue
    #b2df8a a light green
    #33a02c a more vivid darker green
    #fb9a99 pastel pink color
    #e31a1c more vivid pink / red
    #fdbf6f light orange
    #ff7f00 dark orange
    #cab2d6 light purple
    #6a3d9a dark purple
    #ffff99 light yellow
    #b15928 brown

    # OK so I want to use four colors. I think it's nice to have two light ones
    # and two dark ones.

    #maxi(ax)
    tde(axarr, '#33a02c', legend=True)
    sn2003L(axarr, '#1f78b4', legend=True)
    sn1979c(axarr, '#1f78b4', None)
    sn1993J(axarr, '#1f78b4', None)
    sn2011dh(axarr, '#1f78b4', None)
    #grb030329(axarr, '#6a3d9a', legend=True)
    #grb130427A(axarr, '#6a3d9a', None)
    #sn2007bg(axarr, '#1f78b4', None)
    #sn2003bg(axarr, '#1f78b4', None)
    #sn2009bb(axarr, '#e31a1c', legend=True)
    #sn1998bw(axarr, '#e31a1c', None)
    at2018cow(axarr, 'k', None)

    axarr[0].set_ylabel(
            r"Luminosity $\nu L_{\nu}$ [erg\,s$^{-1}$]", 
            fontsize=16)
    axarr[0].set_title("$\\nu > 90\,$GHz", fontsize=14)
    axarr[1].set_title("$\\nu < 10\,$GHz", fontsize=14)
    for ax in axarr:
        ax.tick_params(axis='both', labelsize=14)
        ax.set_xlim(0.3, 2000) 
        ax.set_ylim(1E34, 1E44)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r"Time [days; rest frame]", fontsize=16)
    axarr[1].legend(fontsize=12, loc='upper right')

    plt.show()
    #plt.savefig("lum_evolution.png")
