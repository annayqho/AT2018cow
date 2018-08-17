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


def at2018cow(ax):
    """ Peak of 215.5 GHz light curve """
    d = Planck15.luminosity_distance(z=0.014).cgs.value
    t = 20
    nu = 215.5E9
    lum = nu * 53.32 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    ax.scatter(t, lum, c='k', marker='*', s=400)
    ax.text(t*1.2, lum, 'AT2018cow', 
            verticalalignment='center', fontsize=14)


def tde(ax):
    """  Peak of the 200 GHz light curve: 14.1 mJy, 16.12 days """
    z = 0.354
    d = Planck15.luminosity_distance(z=z).cgs.value
    nu = np.array([4.9E9, 200E9])
    flux = np.array([2.11, 14.1])
    lum = nu * flux * 1e-23 * 1e-3 * 4 * np.pi * d**2
    t = np.array([19.73, 16.12]) / (1+z)
    ax.scatter(
            t[0], lum[0], marker='o', edgecolor='k', facecolor='white', s=100)
    ax.scatter(
            t[1], lum[1], marker='o', edgecolor='k', facecolor='black', s=100)
    ax.plot(
            t, lum, ls='--', c='k')
    ax.text(t[1]*1.3, lum[1], 'SwiftJ1644+57', 
            verticalalignment='center', fontsize=11)


def sn2003L(ax):
    """ Soderberg et al """
    d = 2.8432575937224894e+26
    nu = 8.5E9
    flux = 2.776
    flux_err = 0.051
    lum = nu * flux * 1e-23 * 1e-3 * 4 * np.pi * d**2
    t = 85.4
    ax.scatter(t, lum, marker='s', edgecolor='k', facecolor='white', s=100)
    

def sn1979c(ax):
    """ Weiler 1986 
    This is a IIL """
    d = 5.341805643483106e+25
    nu = 1.4E9
    flux = 10
    lum = nu * flux * 1e-23 * 1e-3 * 4 * np.pi * d**2
    t = 1300 # very approximate
    ax.scatter(t, lum, marker='s', edgecolor='k', facecolor='white', s=100)
    

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
    lum = freq * flux * 1e-23 * 1e-3 * 4 * np.pi * d**2
    # lower freq has no uncertainty
    #lum_err = freq * 4.5 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    ax.scatter(
            t[0], lum[0], marker='s', edgecolor='k', 
            facecolor='white', s=100)
    ax.scatter(t[1], lum[1], marker='s', edgecolor='k', 
            facecolor='black', s=100)
    ax.plot(t, lum, ls='--', c='k')
    ax.text(t[0]/1.5, lum[0], "SN1993J", fontsize=11,
            horizontalalignment='right')


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
    lum = nu * f * 1e-23 * 1e-3 * 4 * np.pi * d**2
    lum_err = nu * f_err * 1e-23 * 1e-3 * 4 * np.pi * d**2
    mcols = ['black', 'white']
    ratios = [1/1.5, 1.5]
    for ii,nuval in enumerate(nu):
        plt.scatter(
                t[ii], lum[ii], marker='s', s=100,
                facecolor=mcols[ii], edgecolor='black')
        plot_limits(ax, t[ii], lum[ii], ratios[ii], 3, mcols[ii])
    plt.plot(t, lum, ls='--', c='k')
    plt.text(
            t[0], lum[0]/2, "SN2011dh", fontsize=11, 
            horizontalalignment='right', verticalalignment='top')


def grb030329(ax):
    """ Sheth et al. 2003 """
    freq = np.array([100E9, 8.5E9])
    z = 0.1686
    t = np.array([3, 14.87]) / (1+z)
    d = Planck15.luminosity_distance(z=z).cgs.value
    f = np.array([70, 19.15])
    lum = freq * f * 1e-23 * 1e-3 * 4 * np.pi * d**2
    f_err = np.array([6, 0.08])
    lum_err = freq * f_err * 1e-23 * 1e-3 * 4 * np.pi * d**2
    mcols = ['black', 'white']
    for ii,nuval in enumerate(freq):
        plt.scatter(
                t[ii], lum[ii], marker='o', s=100, 
                facecolor=mcols[ii], edgecolor='k')
    plt.plot(t, lum, c='k', ls='--')
    plt.text(t[0], lum[0]*2, "GRB030329", fontsize=11,
            horizontalalignment='center')


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
    lum = freq * flux * 1e-23 * 1e-3 * 4 * np.pi * d**2
    flux_err = np.array([0.088, 0.365])
    lum_err = flux_err * freq * 1e-23 * 1e-3 * 4 * np.pi * d**2
    mcols = ['white', 'black']
    for ii,nuval in enumerate(freq):
        plt.scatter(
                t[ii], lum[ii], marker='o', s=100,
                facecolor=mcols[ii], edgecolor='k')
    plot_limits(ax, t[1], lum[1], 1/1.5, 3, mcols[1])
    plt.plot(t, lum, c='k', ls='--')
    plt.text(t[1], lum[1]*3, "GRB130427A", fontsize=11,
            horizontalalignment='center')


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
    lum = nu * 1.490 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    lum_err = nu * 0.104 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    plt.scatter(t, lum, marker='s', s=100, 
            facecolor='white', edgecolor='black')


def sn2003bg(ax):
    """ Soderberg et al.
    Peak is resolved for 22.5, 15, 8.46, 4.86, 1.43
    Again, there are two peaks...
    Let's choose the first peak, 8.46
    """
    nu = 8.46E9
    t = 58
    d = 6.056450393620008e+25
    lum = nu * 51.72 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    lum_err = nu * 1.04 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    plt.scatter(
            t, lum, marker='s', s=100,
            facecolor='white', edgecolor='black')


def sn2009bb(ax):
    """ expl date Mar 19 """
    t = 17
    d = 1.237517263280789e+26
    nu = 8.46E9
    flux = 24.681
    flux_err = 0.066
    lum = nu * flux * 1e-23 * 1e-3 * 4 * np.pi * d**2
    lum_err = nu * flux_err * 1e-23 * 1e-3 * 4 * np.pi * d**2
    plt.scatter(t, lum, marker='o', s=100, edgecolor='k', facecolor='white')



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
    """ SN 1998bw """
    t = 12.4
    d = 1.17E26 # cm
    nu = 150E9
    lum = nu * 39 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    plt.scatter(t, lum, marker='o',
            facecolor='k', edgecolor='k', s=100)

    # SN 2013ak

    # most luminous events are typically SNe IIn (Chevalier 2006)

    # SN Ib/c: Berger 2003, Chevalier & Fransson 2006
    # 2008D


if __name__=="__main__":
    fig, ax = plt.subplots(1, 1, figsize=(8,6), sharex=True, sharey=True)

    at2018cow(ax)
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

    plt.show()
