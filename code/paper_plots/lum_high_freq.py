""" Histogram of luminosities at high frequencies

"High" includes: 215.5, 200 GHz, 99.4 GHz
GRBs
SNe
TDEs
"""


import matplotlib
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import Planck15


def at2018cow():
    """ Peak of 215.5 GHz light curve """
    d = Planck15.luminosity_distance(z=0.014).cgs.value
    t = 20
    nu = 215.5E9
    lum = 53.32 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    return t, nu*lum


def tde():
    """  Peak of the 200 GHz light curve: 14.1 mJy, 16.12 days """
    z = 0.354
    d = Planck15.luminosity_distance(z=z).cgs.value
    nu = 200E9
    lum = 14.1 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    name = 'Swift1644+57'
    t = 16.12 / (1+z)
    #velocity.append(0.5*1.1)
    #density.append(0.2E4)
    return name, t, nu*lum


def sn1993J():
    # This is the peak of the 99.4 GHz light curve
    # There is also a 110 GHz point, but only one,
    # so I can't get the peak.
    t = 66.33
    d = 1.1e25
    freq = 99.4E9
    lum = freq * 22.0 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    lum_err = freq * 4.5 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    return t, lum, lum_err


def sn2011dh():
    # SN 2011dh
    # data at the youngest phase ever of a CC SN: days 3-12
    # M51: d = 8.03 Mpc; expl date May 31.58
    # There are limits... peak was earlier and more luminous
    d = 2.5E25
    nu = 107E9
    t = 4.08
    lum = nu * 4.55 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    lum_err = nu * 0.79 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    return t, lum, lum_err


def grb030329():
    """ Sheth et al. 2003 """
    freq = 100E9
    z = 0.1686
    t = 3 / (1+z)
    d = Planck15.luminosity_distance(z=z).cgs.value
    lum = freq * 70 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    lum_err = freq * 6 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    return t, lum, lum_err


def grb130427A():
    """ Perley et al
    They have data from CARMA/PdBI at 90 GHz (3mm)
    But by the time they caught it, it was fading
    """
    freq = 93E9
    z = 0.340
    t = 0.76900 / (1+z)
    d = Planck15.luminosity_distance(z=z).cgs.value
    lum = freq * 3.416 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    lum_err = 0.365 * freq * 1e-23 * 1e-3 * 4 * np.pi * d**2
    return t, lum, lum_err


def grb():
    """ GRBs and Rel SNe """
    names = []
    t = []
    lum = []

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

    # SN 1998bw
    t.append(12.4)
    d = 1.17E26 # cm
    nu = 150E9
    lum.append(nu * 39 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    names.append("SN 1998bw")
    #velocity.append(1.8)

    # SN 2013ak

    # most luminous events are typically SNe IIn (Chevalier 2006)

    # SN Ib/c: Berger 2003, Chevalier & Fransson 2006
    # 2008D
    
    return names, t, lum 


if __name__=="__main__":
    fig, axarr = plt.subplots(1, 2, figsize=(9,5), sharey=True)

    # freq > 90 GHz
    ax = axarr[0] 
    ax.text(0.1, 0.9, r"$\nu > 90\,$GHz", fontsize=14, transform=ax.transAxes)

    t, lum = at2018cow()
    ax.scatter(t, lum, c='k', marker='*', s=250)
    ax.text(t+2, lum, 'AT2018cow', 
            verticalalignment='center', fontsize=14)

    name, t, lum = tde()
    ax.scatter(
            t, lum, marker='o', edgecolor='k', facecolor='white',
            label="Rel. TDE", s=100)
    ax.text(t+1, lum, 'SwiftJ1644+57', 
            verticalalignment='center', fontsize=11)

    t, lum, lum_err = sn1993J()
    ax.scatter(
        t, lum, marker='o', s = 100, c='k',
        label="Interacting SN")
    ax.text(
            t+4, lum, "SN1993J", fontsize=11, 
            verticalalignment='center')

    t, lum, lum_err = sn2011dh()
    ax.errorbar(
            t, lum, yerr=lum_err, fmt='o', ms=10, c='k')
    ax.annotate('', xy=(t/1.5,lum), xytext=(t,lum),
            arrowprops=dict(
                facecolor='black', headwidth=10, width=1, headlength=7))
    ax.annotate('', xy=(t,1E38), xytext=(t,lum),
            arrowprops=dict(
                facecolor='black', headwidth=10, width=1, headlength=7))
    ax.text(
            t+0.3, lum, "SN2011dh", fontsize=11,
            horizontalalignment='left', verticalalignment='center')


    # GRBs
    t, lum, lum_err = grb030329()
    ax.errorbar(t, lum, yerr=lum_err, fmt='s', ms=10,
                 mec='k', mfc='white', label="GRB")
    ax.text(
            t*1.1, lum, "GRB030329", fontsize=11,
            horizontalalignment='left', verticalalignment='center')


    t, lum, lum_err = grb130427A()
    ax.text(
            t+0.06, lum, "GRB130427A", fontsize=11,
            horizontalalignment='left', verticalalignment='center')
    ax.annotate('', xy=(t/1.5,lum), xytext=(t,lum),
            arrowprops=dict(
                facecolor='white', headwidth=10, width=1, headlength=7))
    ax.annotate('', xy=(t,3E42), xytext=(t,lum),
            arrowprops=dict(
                facecolor='white', headwidth=10, width=1, headlength=7))
    ax.errorbar(t, lum, yerr=lum_err, fmt='s', ms=10,
                 mec='k', mfc='white')


    # Make pretty plot
    ax.set_ylabel(
            r"Peak Luminosity $\nu L_{\nu}$ [erg\,s$^{-1}$]", 
            fontsize=16)
    ax.set_xlabel(r"Time to Peak [days; rest frame]", fontsize=16)

    #ax.set_ylabel("Number of Objects", fontsize=16)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.3, 230) 
    #ax.get_yaxis().set_visible(False)
    ax.set_ylim(1E37, 1E44)
    ax.legend(fontsize=12, loc='center left')
    ax.tick_params(axis='both', labelsize=14)

    plt.tight_layout()
    plt.show()
    #plt.savefig("early_nu_lnu.png")
