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
    d = Planck15.luminosity_distance(z=0.354).cgs.value
    nu = 200E9
    lum = 14.1 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    name = 'Swift1644+57'
    t = 16.12
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


def grb():
    """ GRBs and Rel SNe """
    names = []
    t = []
    lum = []

    names.append("GRB 030329")
    freq = 100E9
    t.append(12)
    d = Planck15.luminosity_distance(z=0.1686).cgs.value
    lum.append(freq * 25.8 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    # 2E17 at 15 days (Berger et al) corresponds to
    #V = 2E17 / (15*86400)
    #velocity.append(5.14)
    #density.append(1.8)

    names.append("GRB 130427A")
    freq = 100E9
    t.append(10.35)
    d = Planck15.luminosity_distance(z=0.340).cgs.value
    lum.append(freq * 0.368 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    #velocity.append(130)
    #density.append()

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
    # you could populate it with GRBs, some supernovae, and TDEs. 
    # I would do two symbols: one for measurements at >80 GHz, 
    # and one for all others. 

    fig, ax = plt.subplots(1, 1, figsize=(7,5))
    t, lum = at2018cow()
    plt.scatter(t, lum, c='k', marker='*', s=200)
    plt.text(t, lum, 'AT2018cow', verticalalignment='bottom', fontsize=14)

    name, t, lum = tde()
    plt.scatter(
            t, lum, marker='o', edgecolor='k', facecolor='white',
            label="Rel. TDE", s=100)

    t, lum, lum_err = sn1993J()
    plt.scatter(
        t, lum, marker='o', s = 100, c='k',
        label="Interacting SN")
    plt.text(t, lum, "SN1993J")

    t, lum, lum_err = sn2011dh()
    plt.errorbar(
            t, lum, yerr=lum_err, fmt='o', ms=10, c='k')
    plt.annotate('', xy=(0.5,0.5), xytext=(0.05, 0.05),
            #xycoords='data',
            textcoords='figure fraction')
            #arrowprops=dict(facecolor='black', shrink=0.005))
    plt.text(t, lum, "SN2011dh")

    names, t, lum = grb()
    for ii,name in enumerate(names):
        if ii == 0:
            plt.scatter(t[ii], lum[ii], marker='s', s=100,
                edgecolor='k', facecolor='white', label="GRB")
            #plt.text(t[ii], lum[ii], name)
        else:
            plt.scatter(t[ii], lum[ii], marker='s', s=100,
                edgecolor='k', facecolor='white')
            #plt.text(t[ii], lum[ii], name)

    ax.set_ylabel(
            r"Peak Luminosity $\nu L_{\nu}$ [erg\,s$^{-1}$]", 
            fontsize=16)
    ax.set_xlabel(r"Time to Peak [days]", fontsize=16)

    #ax.set_ylabel("Number of Objects", fontsize=16)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(3, 100) 
    #ax.get_yaxis().set_visible(False)
    ax.set_ylim(1E37, 1E44)
    plt.tight_layout()
    plt.legend(fontsize=12)
    ax.tick_params(axis='both', labelsize=14)
    plt.show()
    #plt.savefig("early_nu_lnu.png")
