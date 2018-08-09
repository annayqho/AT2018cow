""" Histogram of luminosities at high frequencies

Vikram suggests: peak measured radio luminosity (F_nu) on the y-axis
Time of measurement on the x-axis

GRBs
SNe
TDEs
Two symbols, one for measurements at > 80 GHz, one for all others
80 GHz symbols should be brighter
Log x-axis limited to 50 days

Procedure:
    1. For each event, take data between 10 and 20 days
    2. Choose points between 93 and 215.5 GHz
    3. Measure luminosity
    4. Histogram of luminosities
"""


import matplotlib
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import Planck15
from plot_spec import get_data


def get_vals():
    lum = []
    velocity = []
    names = []
    density = []

    ## AT2018cow

    # On Day 10, peak is 144 GHz at 92 mJy
    d = Planck15.luminosity_distance(z=0.014).cgs.value

    # 230 GHz flux varies from 15-45 mJy,
    # but that's above 200 GHz

    # Let's say 215.5. Peak of that 
    # in this date range is 53.32
    # that's day 20, which is very close to day 22

    lum.append(215.5e9 * 53.32 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    names.append('AT2018cow')
    velocity.append(0.17*1.01)
    density.append(5.3E4)

    ## swift1644
    # On Day 10, peak is 40 mJy at 250 GHz
    # On Day 15, peak is 30 mJy at 140 GHz
    # actual measurement: 16.12 days, 200 GHz, 14.1
    d = Planck15.luminosity_distance(z=0.354).cgs.value
    lum.append(200e9 * 14.1 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    names.append('Swift1644+57')
    # using the values measured on Day 15
    velocity.append(0.5*1.1)
    density.append(0.2E4)


    ## SN 1993J
    # 19 days
    d = 1.1e25
    lum.append(110e9 * 10.0 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    names.append('SN 1993J') 
    # velocity
    R = d * 74.6 * 1E-6 / 206265
    V = R/(19*86400)
    velocity.append(0.08)


    # GRB 030329
    # 100 GHz, Sheth et al. 2003, 25.8 on Day 12
    d = Planck15.luminosity_distance(z=0.1686).cgs.value
    lum.append(100e9 * 25.8 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    names.append("GRB 030329")
    # 2E17 at 15 days (Berger et al) corresponds to
    V = 2E17 / (15*86400)
    velocity.append(5.14)
    density.append(1.8)

    # GRB 130427A
    # 100 GHz, Perley et al. 0.368 mJy
    d = Planck15.luminosity_distance(z=0.340).cgs.value
    lum.append(100e9 * 0.368 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    names.append("GRB 130427A")
    velocity.append(130)
    #density.append()

    # GRB 070125, GRB 100418A, GRB 970508, 090423
    # de Ugarte Postigo et al 2012
    # 070125: Castro-Tirado et al
    d = Planck15.luminosity_distance(z=1.55).cgs.value
    lum.append(100e9 * 1.23 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    names.append("GRB 070125")
    velocity.append(100) # placeholder
    
    # GRB 100418A
    d = Planck15.luminosity_distance(z=0.62).cgs.value
    lum.append(106e9 * 1.13 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    names.append("GRB 100418A")
    velocity.append(300) # temporary placeholder

    # SNe II: Markarian 297 (Yin & Heeschen 1991)
    # FIRST J121550.2 + 130654 (Levinson et al. 2002; Gal-Yam et al. 2006)
    # SN 2008iz (Brunthaler et al. 2009, 2010)

    # SN 2011dh
    # data at the youngest phase ever of a CC SN: days 3-12
    # M51: d = 8.03 Mpc; expl date May 31.58
    d = 2.5E25
    lum.append(93e9 * 1.61 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    names.append("SN 2011dh")
    velocity.append(0.15)
    density.append(5E3)

    # SN 1998bw
    d = 1.17E26 # cm
    lum.append(150E9 * 39 * 1e-23 * 1e-3 * 4 * np.pi * d**2)
    names.append("SN 1998bw")
    velocity.append(1.8)

    # SN 2013ak

    # most luminous events are typically SNe IIn (Chevalier 2006)

    # SN Ib/c: Berger 2003, Chevalier & Fransson 2006
    # 2008D
    
    return names, lum, velocity, density



if __name__=="__main__":
    fig, ax = plt.subplots(1, 1, figsize=(5,5))
    names, lum, v, n = get_vals()

    plt.scatter(v, lum)
    for ii,val in enumerate(names):
        plt.text(v[ii], lum[ii], val, verticalalignment='bottom')

    #for ii,val in enumerate(lum):
    #    plt.axvline(x=val, c='k', lw=8.0)
    #    plt.text(val*1.05, ii*0.1, names[ii], fontsize=14)

    ax.set_ylabel("$\\nu L_\\nu$ [erg\,s$^{-1}$]", fontsize=16)
    ax.set_xlabel("$\\beta \\Gamma$")

    #ax.set_ylabel("Number of Objects", fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.set_xlim(1E37, 1E44) 
    #ax.get_yaxis().set_visible(False)
    #ax.set_ylim(0.5, 1000)
    plt.tight_layout()
    plt.show()
    #plt.savefig("early_nu_lnu.png")
