""" Histogram of luminosities at high frequencies

GRBs
SNe
TDEs
Two symbols, one for measurements at > 80 GHz, one for all others
80 GHz symbols should be brighter
Log x-axis limited to 50 days
"""


import matplotlib
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import Planck15
from plot_spec import get_data


def at2018cow():
    """ Peak of 215.5 GHz light curve """
    d = Planck15.luminosity_distance(z=0.014).cgs.value
    t = 20
    nu = 215.5
    lum = 53.32 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    # v = 0.17 * 1.01
    # density = 5.3E4
    return t, lum


def tde():
    """  Peak of the 200 GHz light curve: 14.1 mJy, 16.12 days """
    d = Planck15.luminosity_distance(z=0.354).cgs.value
    nu = 200E9
    lum = 14.1 * 1e-23 * 1e-3 * 4 * np.pi * d**2
    name = 'Swift1644+57'
    t = 16.12
    #velocity.append(0.5*1.1)
    #density.append(0.2E4)
    return name, t, lum


def sn():
    """ Non-relativistic supernovae """
    names = []
    t = []
    lum = []

    names.append("SN1993J") 
    t.append(19)
    d = 1.1e25
    freq = 110E9
    lum.append(10.0 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    # velocity
    R = d * 74.6 * 1E-6 / 206265
    V = R/(19*86400)
    #velocity.append(0.08)


    # SN 2011dh
    # data at the youngest phase ever of a CC SN: days 3-12
    # M51: d = 8.03 Mpc; expl date May 31.58
    names.append("SN2011dh")
    d = 2.5E25
    # nu = 93E9
    t.append(10.95)
    lum.append(1.61 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    #velocity.append(0.15)
    #density.append(5E3)


    return names, t, lum


def grb():
    """ GRBs and Rel SNe """
    names = []
    t = []
    lum = []

    names.append("GRB 030329")
    freq = 100
    t.append(12)
    d = Planck15.luminosity_distance(z=0.1686).cgs.value
    lum.append(25.8 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    # 2E17 at 15 days (Berger et al) corresponds to
    #V = 2E17 / (15*86400)
    #velocity.append(5.14)
    #density.append(1.8)

    names.append("GRB 130427A")
    freq = 100
    t.append(10.35)
    d = Planck15.luminosity_distance(z=0.340).cgs.value
    lum.append(0.368 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    #velocity.append(130)
    #density.append()

    # first GRB detected by ALMA
    #names.append("GRB 110715A")
    # freq = 345 GHz

    # recent LLGRB
    names.append("GRB 171205A")
    freq = 92
    t.append(6)
    d = Planck15.luminosity_distance(z=0.037).cgs.value
    lum.append(28 * 1e-23 * 1e-3 * 4 * np.pi * d**2)

    # GRB 070125, GRB 100418A, GRB 970508, 090423
    # de Ugarte Postigo et al 2012
    # 070125: Castro-Tirado et al
    names.append("GRB 070125")
    d = Planck15.luminosity_distance(z=1.55).cgs.value
    t.append(14.79)
    nu = 100E9
    lum.append(1.23 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    
    names.append("GRB 100418A")
    t.append(12.328)
    # nu = 106E9
    d = Planck15.luminosity_distance(z=0.62).cgs.value
    lum.append(1.13 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    #velocity.append(300) # temporary placeholder

    # SN 1998bw
    t.append(12.4)
    d = 1.17E26 # cm
    # nu = 150E9
    lum.append(39 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
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

    names, t, lum = sn()
    for ii,name in enumerate(names):
        if ii == 0:
            plt.scatter(
                t[ii], lum[ii], marker='o', s = 100, c='k',
                label="Interacting SN")
        else:
            plt.scatter(t[ii], lum[ii], marker='o', s=100, c='k')

    names, t, lum = grb()
    for ii,name in enumerate(names):
        if ii == 0:
            plt.scatter(t[ii], lum[ii], marker='s', s=100,
                edgecolor='k', facecolor='white', label="GRB")
        else:
            plt.scatter(t[ii], lum[ii], marker='s', s=100,
                edgecolor='k', facecolor='white')

    ax.set_ylabel(r"$L_{\nu}$ [erg\,s$^{-1}$\,Hz$^{-1}$]", fontsize=16)
    ax.set_xlabel(r"$\Delta t$ [days]", fontsize=16)

    #ax.set_ylabel("Number of Objects", fontsize=16)
    #ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(5, 25) 
    #ax.get_yaxis().set_visible(False)
    ax.set_ylim(1E25, 1E33)
    plt.tight_layout()
    plt.legend(fontsize=12)
    ax.tick_params(axis='both', labelsize=14)
    plt.show()
    #plt.savefig("early_nu_lnu.png")
