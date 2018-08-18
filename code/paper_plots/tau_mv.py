""" 
A Tau Mv plot showing that we expect there to have been
more luminous millimeter transients
"""


import matplotlib
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import Planck15


def plot_limits(ax, x, y, ratiox, ratioy, col):
    """ Plot two arrows from the point """
    ax.annotate('', xy=(x*ratiox, y), xytext=(x, y),
            arrowprops=dict(
                facecolor=col, headwidth=10, width=1, headlength=7))
    ax.annotate('', xy=(x, y*ratioy), xytext=(x, y),
            arrowprops=dict(
                facecolor=col, headwidth=10, width=1, headlength=7))


def plot_point(ax, d, nu, t, f, marker, name=None):
    """ Plot a single point
    If nu > 90 GHz, make it black
    If nu < 10 GHz, make it white
    
    Parameters
    ----------
    f: flux in mJy
    """
    if nu < 10E9:
        col = 'white'
    elif nu > 90E9:
        col = 'black'
    else:
        print("unknown freq range")
    if marker=='*':
        s=400
    else:
        s=100
    lum = nu * f * 1e-23 * 1e-3 * 4 * np.pi * d**2
    ax.scatter(t, lum, facecolor=col, edgecolor='k', marker=marker, s=s)
    if name:
        if name=='AT2018cow':
            fs = 14
        else:
            fs = 11
        ax.text(t*1.2, lum, name, 
                verticalalignment='center', fontsize=fs)
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
    """ Peak of 215.5 GHz light curve """
    d = Planck15.luminosity_distance(z=0.014).cgs.value
    t = 20
    nu = 215.5E9
    plot_point(ax, d, nu, t, 53.32, '*', name='AT2018cow')


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
    """  Peak of the 200 GHz light curve: 14.1 mJy, 16.12 days """
    z = 0.354
    d = Planck15.luminosity_distance(z=z).cgs.value
    nu = np.array([4.9E9, 200E9])
    flux = np.array([2.11, 14.1])
    t = np.array([19.73, 16.12]) / (1+z)
    plot_points(ax, d, nu, t, flux, marker='o', name='SwiftJ1644+57')


def sn2003L(ax):
    """ Soderberg et al """
    d = 2.8432575937224894e+26
    nu = 8.5E9
    flux = 2.776
    flux_err = 0.051
    t = 85.4
    plot_point(ax, d, nu, t, flux, 's')
    

def sn1979c(ax):
    """ Weiler 1986 
    This is a IIL """
    d = 5.341805643483106e+25
    nu = 1.4E9
    flux = 10
    t = 1300 # very approximate
    plot_point(ax, d, nu, t, flux, 's')
    

def sn1993J(ax):
    """ SN 1993J from Weiler et al. 
    This is the peak of the 99.4 GHz light curve
    There is also a 110 GHz point, but only one,
    so I can't get the peak.
    1.4 GHz / 20 cm, 4.9 GHz / 6 cm
    values come from the best-fit model,
    but by eye they are clearly pretty close
    """
    t = np.array([133, 66.33])
    d = 1.1e25
    freq = np.array([4.9E9, 99.4E9])
    flux = np.array([96.9, 22])
    # lower freq has no uncertainty
    #lum_err = freq * 4.5 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    plot_points(ax, d, freq, t, flux, 's', 'SN1993J')


def sn2011dh(ax):
    """ SN 2011dh
    M51: d = 8.03 Mpc; expl date May 31.58
    The high-freq points are limits. The peak was earlier and more luminous.
    The low-freq points are limits. The peak was later and more luminous.
    """
    d = 2.5E25
    nu = np.array([107E9, 8.5E9])
    t = np.array([4.08, 11.98])
    f = np.array([4.55, 3.15])
    f_err = np.array([0.79, 0.218])
    lums = plot_points(ax, d, nu, t, f, 's', 'SN2011dh')
    ratios = [1/1.5, 1.5]
    for ii,nuval in enumerate(nu):
        plot_limits(ax, t[ii], lums[ii], ratios[ii], 3, 'k')


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
    fig, ax = plt.subplots(1, 1, figsize=(8,6), sharex=True, sharey=True)

    at2018cow(ax)
    maxi(ax)
    tde(ax)
    sn2003L(ax)
    sn1979c(ax)
    sn1993J(ax)
    sn2011dh(ax)
    grb030329(ax)
    grb130427A(ax)
    sn2007bg(ax)
    sn2003bg(ax)
    sn2009bb(ax)
    sn1998bw(ax)

    ax.tick_params(axis='both', labelsize=14)
    #ax.set_xlim(0.3, 2000) 
    #ax.set_ylim(1E36, 1E44)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(
            r"Peak Luminosity $\nu L_{\nu}$ [erg\,s$^{-1}$]", 
            fontsize=16)
    ax.set_xlabel(r"Time to Peak [days; rest frame]", fontsize=16)
    #ax.legend(fontsize=12, loc='center left')

    #plt.show()
    plt.savefig("mm_tau_lum.png")
