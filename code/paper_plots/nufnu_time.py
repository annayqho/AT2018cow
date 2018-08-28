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
from astropy.cosmology import Planck15
from get_radio import *


def plot_limits(ax, x, y, ratiox, ratioy, col):
    """ Plot two arrows from the point """
    ax.annotate('', xy=(x*ratiox, y), xytext=(x, y),
            arrowprops=dict(
                facecolor=col, headwidth=10, width=1, headlength=7))
    ax.annotate('', xy=(x, y*ratioy), xytext=(x, y),
            arrowprops=dict(
                facecolor=col, headwidth=10, width=1, headlength=7))


def plot_line(ax, d, nu, t, f, name, label):
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
    if name=='AT2018cow':
        marker='*'
        s=150
        col='k'
    else:
        if label=='SN':
            marker='o'
            col='k'
            s=50
        elif label=='GRB':
            marker='o'
            col='white'
            s=50
        elif label=='Rel. SN':
            marker='s'
            col='k'
            s=50
        elif label=='Rel. TDE':
            marker='s'
            col='white'
            s=50
    ax.scatter(
            t, lum, facecolor=col, edgecolor='black', 
            marker=marker, s=s)
    ax.plot(t, lum, c='black', ls='-')
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


def at2018cow(ax):
    """ 231.5 GHz light curve and 9 GHz light curve """
    tel, freq, days, flux, eflux_form, eflux_sys = get_data_all()
    d = Planck15.luminosity_distance(z=0.014).cgs.value
    choose = np.logical_and(tel=='SMA', freq==231.5)
    nu = 231.5E9
    lum = plot_line(
            ax[0], d, nu, days[choose], flux[choose], 'AT2018cow', None)

    ax[0].text(days[-1]*1.05, lum[-1], 'AT2018cow', fontsize=11,
            verticalalignment='center',
            horizontalalignment='left')
    choose = np.logical_and(tel=='ATCA', freq==9)
    nu = 9E9
    plot_line(
            ax[1], d, nu, days[choose], flux[choose], 'AT2018cow', None)



def maxi(ax):
    """ data from Dillon """
    d = 4.62E26
    t = 4*365
    nu = 6E9
    f = 6.16
    lum = plot_point(ax, d, nu, t, f, '*')
    plt.text(t/1.2, lum, 'MAXI 140814A', fontsize=11,
            horizontalalignment='right', verticalalignment='center')



def tde(ax):
    """  Plot the 225 GHz light curve from the SMA
    
    Plot the 4.9 GHz light curve from the VLA
    """
    z = 0.354
    d = Planck15.luminosity_distance(z=z).cgs.value
    nu = 225E9
    t = np.array([7.20, 7.86, 14.10, 15.12, 17.11, 18.13]) / (1+z)
    flux = np.array([14.9, 11.7, 13.3, 9.9, 8.2, 8.3])
    plot_line(ax[0], d, nu, t, flux, 'SwiftJ1644+57', 'Rel. TDE')
    nu = 4.9E9
    t = np.array([0.82, 1.71, 1.95, 2.74, 3.73, 4.72, 6.74, 11.93, 19.73]) / (1+z)
    flux = np.array([0.25, 0.34, 0.34, 0.61, 0.82, 1.48, 1.47, 1.80, 2.11])
    plot_line(ax[1], d, nu, t, flux, 'SwiftJ1644+57', 'Rel. TDE')


def sn2003L(ax):
    """ Soderberg et al
    Values at 8.5 GHz """
    d = 2.8432575937224894e+26
    nu = 8.5E9
    t = np.array([25.2, 27.3, 28.3, 29.4, 30.3, 31.2, 32.3, 
        36.3, 38.3, 41.6, 44.5, 46.4, 48.4, 
        53.3, 55.3, 60.3])
    flux = np.array([743, 810, 848, 966, 1051, 947, 883, 1147, 
        1211, 1483, 1448, 1413, 1546, 1854, 1969, 2283])
    plot_line(ax[1], d, nu, t, flux, 's', 'SN2003L', 'left')
    

def sn1979c(ax):
    """ Weiler 1986 
    This is a IIL """
    d = 5.341805643483106e+25
    nu = 1.4E9
    t = np.array([437, 594, 631, 727, 822, 914, 1026, 1156, 1212])
    flux = np.array([0.2, 2.1, 2.5, 4.4, 7.1, 9.8, 10.2, 12.2, 10.2])
    plot_line(ax[1], d, nu, t, flux, 's', 'SN1979c', 'left')
    

def sn1993J(ax):
    """ SN 1993J from Weiler et al. 
    This is the peak of the 99.4 GHz light curve
    There is also a 110 GHz point, but only one,
    so I can't get the peak.
    1.4 GHz / 20 cm, 4.9 GHz / 6 cm
    values come from the best-fit model,
    but by eye they are clearly pretty close
    """
    d = 1.1e25
    freq = 4.9E9
    t = np.array([16.53, 17.07, 19.02, 22.25, 22.97, 24.51, 25.48, 26.88, 27.87, 28.99, 29.85, 35.71, 40.11, 53.17, 75.06, 102.76])
    flux = np.array([0.327, 0.280, 0.360, 0.880, 0.870, 1.290, 1.700, 1.930, 2.050, 2.640, 3.260, 6, 10.510, 25.819, 56.330, 102.670])
    plot_line(ax[1], d, freq, t, flux, 's', 'SN1993J', 'left')
    freq = 99.4E9
    t = np.array([14.39, 17.39, 24.38, 33.28, 43.39, 63.32, 66.33, 82.23, 97.17, 168.80, 195.61, 231.67])
    flux = np.array([18, 17, 20, 17, 23, 19, 22, 14, 16, 13, 8, 8])
    plot_line(ax[0], d, freq, t, flux, 's', 'SN1993J', 'left')


def sn2011dh(ax):
    """ SN 2011dh
    M51: d = 8.03 Mpc; expl date May 31.58
    The high-freq points are limits. The peak was earlier and more luminous.
    The low-freq points are limits. The peak was later and more luminous.
    """
    d = 2.5E25
    nu = 107E9
    t = np.array([4.08, 5.22])
    flux = np.array([4.55, 3.66])
    plot_line(ax[0], d, nu, t, flux, 's', 'SN2011dh', 'right')
    nu = 8.5E9
    t = np.array([5.01, 7.12, 9.02, 11.98])
    flux = np.array([0.455, 1.06, 1.58, 3.15])
    plot_line(ax[1], d, nu, t, flux, 's', 'SN2011dh', 'left')


def grb030329(ax):
    """ Sheth et al. 2003 """
    freq = np.array([100E9, 8.5E9])
    z = 0.1686
    t = np.array([3, 14.87]) / (1+z)
    d = Planck15.luminosity_distance(z=z).cgs.value
    f = np.array([70, 19.15])
    f_err = np.array([6, 0.08])
    plot_points(ax, d, freq, t, f, 'o', 'GRB030329')


def grb130427A(ax):
    """ Perley et al
    They have data from CARMA/PdBI at 90 GHz (3mm)
    But by the time they caught it, it was fading
    """
    freq = np.array([5.10E9, 93E9])
    z = 0.340
    t = np.array([2.03854, 0.76900]) / (1+z)
    d = Planck15.luminosity_distance(z=z).cgs.value
    flux = np.array([1.760, 3.416])
    flux_err = np.array([0.088, 0.365])
    lums = plot_points(ax, d, freq, t, flux, 'o', 'GRB130427A')
    plot_limits(ax, t[1], lums[1], 1/1.4, 2.5, 'k')


def sn2007bg(ax):
    """ Salas et al. 2013
    Peak is resolved for 4.86, 8.46 GHz
    Let's choose the one where it's more clearly resolved
    It's a bit weird because there are two peaks, but let's
    choose the first one because it neglects subsequent interaction
    """
    nu = 8.46E9
    t = 55.9
    d = Planck15.luminosity_distance(z=0.0346).cgs.value
    f = 1.490
    f_err = 0.104
    plot_point(ax, d, nu, t, f, 's')


def sn2003bg(ax):
    """ Soderberg et al.
    Peak is resolved for 22.5, 15, 8.46, 4.86, 1.43
    Again, there are two peaks...
    Let's choose the first peak, 8.46
    """
    nu = 8.46E9
    t = 58
    d = 6.056450393620008e+25
    f = 51.72
    f_err = 1.04
    plot_point(ax, d, nu, t, f, 's')


def sn2009bb(ax):
    """ expl date Mar 19 """
    t = 17
    d = 1.237517263280789e+26
    nu = 8.46E9
    flux = 24.681
    flux_err = 0.066
    plot_point(ax, d, nu, t, flux, 'o')


def grb():
    """ GRBs and Rel SNe """
    # first GRB detected by ALMA
    #names.append("GRB 110715A")
    # freq = 345 GHz

    # recent LLGRB
    names.append("GRB 171205A")
    freq = 92E9
    t.append(6)
    d = Planck15.luminosity_distance(z=0.037).cgs.value
    lum.append(freq * 28 * 1e-23 * 1e-3 * 4 * np.pi * d**2)

    # GRB 070125, GRB 100418A, GRB 970508, 090423
    # de Ugarte Postigo et al 2012
    # 070125: Castro-Tirado et al
    names.append("GRB 070125")
    d = Planck15.luminosity_distance(z=1.55).cgs.value
    t.append(14.79)
    nu = 100E9
    lum.append(nu * 1.23 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    
    names.append("GRB 100418A")
    t.append(12.328)
    nu = 106E9
    d = Planck15.luminosity_distance(z=0.62).cgs.value
    lum.append(nu * 1.13 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    #velocity.append(300) # temporary placeholder


def sn1998bw(ax):
    """ SN 1998bw
    
    This is a bit complicated because there are two peaks in the light curve,
    but I am using a frequency that has a main peak rather than a frequecy
    with two clear distinct peaks...
    """
    t = np.array([12.4, 33])
    d = 1.17E26 # cm
    nu = np.array([150E9, 2.3E9])
    f = np.array([39, 34])
    lums = plot_points(ax, d, nu, t, f, 'o', 'SN1998bw')
    plot_limits(ax, t[0], lums[0], 0, 2, 'k')
    plot_limits(ax, t[0], lums[0], 1/1.3, 0, 'k')
    plot_limits(ax, t[0], lums[0], 1.3, 0, 'k')


def othersn(ax):
    """ 
    SN 2013ak
    most luminous events are typically SNe IIn (Chevalier 2006)
    SN Ib/c: Berger 2003, Chevalier & Fransson 2006
    2008D
    """


if __name__=="__main__":
    fig, axarr = plt.subplots(1, 2, figsize=(10,5), sharex=True)
    props = dict(boxstyle='round', facecolor='white')

    at2018cow(axarr)
    #maxi(ax)
    tde(axarr)
    #sn2003L(axarr)
    #sn1979c(axarr)
    #sn1993J(axarr)
    #sn2011dh(axarr)
    # grb030329(ax)
    # grb130427A(ax)
    # sn2007bg(ax)
    # sn2003bg(ax)
    # sn2009bb(ax)
    # sn1998bw(ax)

    axarr[0].set_ylabel(
            r"Luminosity $\nu L_{\nu}$ [erg\,s$^{-1}$]", 
            fontsize=16)
    axarr[0].text(
            0.95, 0.9, "$\\nu > 90\,$GHz", 
            transform=axarr[0].transAxes, fontsize=14, bbox=props,
            horizontalalignment='right')
    axarr[1].text(
            0.95, 0.9, "$\\nu < 10\,$GHz", 
            transform=axarr[1].transAxes, fontsize=14, bbox=props,
            horizontalalignment='right')
    for ax in axarr:
        ax.tick_params(axis='both', labelsize=14)
        ax.set_xlim(0.3, 2000) 
        #ax.set_ylim(1E36, 1E44)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r"Time [days; rest frame]", fontsize=16)
        #ax.legend(fontsize=12, loc='center left')

    plt.show()
    #plt.savefig("lum_evolution.png")
