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
    nu = np.array([4.9E9, 200E9])
    flux = np.array([2.11, 14.1])
    lum = nu * flux * 1e-23 * 1e-3 * 4 * np.pi * d**2
    t = np.array([19.73, 16.12]) / (1+z)
    return t, lum


def sn2003L():
    """ Soderberg et al """
    d = 2.8432575937224894e+26
    nu = 8.5E9
    flux =2.776
    flux_err = 0.051
    lum = nu * flux * 1e-23 * 1e-3 * 4 * np.pi * d**2
    t = 85.4
    return t, lum
    

def sn1979c():
    """ Weiler 1986 
    This is a IIL """
    d = 5.341805643483106e+25
    nu = 1.4E9
    flux = 10
    lum = nu * flux * 1e-23 * 1e-3 * 4 * np.pi * d**2
    t = 1300 # very approximate
    return t, lum
    


def sn1993J():
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
    return t, lum


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
    freq = np.array([5.10E9, 93E9])
    z = 0.340
    t = np.array([2.03854, 0.76900]) / (1+z)
    d = Planck15.luminosity_distance(z=z).cgs.value
    flux = np.array([1.760, 3.416])
    lum = freq * flux * 1e-23 * 1e-3 * 4 * np.pi * d**2
    flux_err = np.array([0.088, 0.365])
    lum_err = flux_err * freq * 1e-23 * 1e-3 * 4 * np.pi * d**2
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


def sn2007bg():
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
    return t, lum, lum_err


def sn2003bg():
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
    return t, lum, lum_err


def sn2009bb():
    """ expl date Mar 19 """
    t = 17
    d = 1.237517263280789e+26
    nu = 8.46E9
    flux = 24.681
    flux_err = 0.066
    lum = nu * flux * 1e-23 * 1e-3 * 4 * np.pi * d**2
    lum_err = nu * flux_err * 1e-23 * 1e-3 * 4 * np.pi * d**2
    return t, lum, lum_err


def hifreq(ax):
    """ 
    Plot the high-frequency part 

    Parameters
    ----------
    ax: the axis to make this plot on
    """
    ax.text(0.1, 0.9, r"$\nu > 90\,$GHz", fontsize=14, transform=ax.transAxes)

    t, lum = at2018cow()
    ax.scatter(t, lum, c='k', marker='*', s=250)
    ax.text(t+2, lum, 'AT2018cow', 
            verticalalignment='center', fontsize=14)

    t, lum = tde()
    ax.scatter(
            t[1], lum[1], marker='o', edgecolor='k', facecolor='white',
            label="Rel. TDE", s=100)
    ax.text(t[1]*1.3, lum[1], 'SwiftJ1644+57', 
            verticalalignment='center', fontsize=11)

    t, lum = sn1993J()
    ax.scatter(
        t[1], lum[1], marker='o', s = 100, c='k',
        label="Interacting SN")
    ax.text(
            t[1]+4, lum[1], "SN1993J", fontsize=11, 
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
    ax.errorbar(t[1], lum[1], yerr=lum_err[1], fmt='s', ms=10,
                 mec='k', mfc='white')
    ax.text(
            t[1]+0.06, lum[1], "GRB130427A", fontsize=11,
            horizontalalignment='left', verticalalignment='center')
    ax.annotate('', xy=(t[1]/1.5,lum[1]), xytext=(t[1],lum[1]),
            arrowprops=dict(
                facecolor='white', headwidth=10, width=1, headlength=7))
    ax.annotate('', xy=(t[1],3E42), xytext=(t[1],lum[1]),
            arrowprops=dict(
                facecolor='white', headwidth=10, width=1, headlength=7))


    # Make pretty plot
    ax.set_ylabel(
            r"Peak Luminosity $\nu L_{\nu}$ [erg\,s$^{-1}$]", 
            fontsize=16)
    ax.set_xlabel(r"Time to Peak [days; rest frame]", fontsize=16)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.3, 230) 
    #ax.legend(fontsize=12, loc='center left')
    ax.tick_params(axis='both', labelsize=14)



if __name__=="__main__":
    fig, axarr = plt.subplots(1, 2, figsize=(9,5), sharex=True, sharey=True)

    ### HIGH-FREQUENCY PLOT
    ax = axarr[0] 
    hifreq(ax)

    ### LOW-FREQUENCY PLOT
    ax = axarr[1]
    ax.text(0.1, 0.9, r"$\nu < 10\,$GHz", fontsize=14, transform=ax.transAxes)
    t, lum, lum_err = sn2007bg()
    ax.scatter(
        t, lum, marker='o', s = 100, c='k')
    # ax.text(
    #         t*1.2, lum, "SN2007bg", fontsize=11, 
    #         verticalalignment='center')

    t, lum, lum_err = sn2003bg()
    ax.scatter(
        t, lum, marker='o', s = 100, c='k')
    # ax.text(
    #         t*1.2, lum, "SN2003bg", fontsize=11, 
    #         verticalalignment='center')

    t, lum, lum_err = grb130427A()
    # lo first, hi second
    ax.errorbar(t[0], lum[0], yerr=lum_err[0], fmt='s', ms=10,
                 mec='k', mfc='white')
    ax.text(
            t[0]*1.2, lum[0], "GRB130427A", fontsize=11, 
            verticalalignment='center')

    t, lum, lum_err = sn2009bb()
    ax.errorbar(t, lum, yerr=lum_err, fmt='s', ms=10,
                 mec='k', mfc='white')
    # ax.text(
    #         t, lum/1.5, "SN2009bb", fontsize=11, horizontalalignment='center',
    #         verticalalignment='top')
    ax.annotate('', xy=(t/1.5,lum), xytext=(t,lum),
            arrowprops=dict(
                facecolor='white', headwidth=10, width=1, headlength=7))
    ax.annotate('', xy=(t,lum*3), xytext=(t,lum),
            arrowprops=dict(
                facecolor='white', headwidth=10, width=1, headlength=7))


    t, lum = sn1993J()
    print(t[0], lum[0])
    ax.scatter(
        t[0], lum[0], marker='o', s = 100, c='k')
    ax.text(
            t[0]/1.2, lum[0], "SN1993J", fontsize=11, 
            verticalalignment='center',
            horizontalalignment='right')
    
    t, lum = sn2003L()
    ax.scatter(
        t, lum, marker='o', s = 100, c='k')
    #ax.text(
    #        t/1.2, lum, "SN2003L", fontsize=11, 
    #        verticalalignment='center',
    #        horizontalalignment='right')


    t, lum = sn1979c()
    ax.scatter(
        t, lum, marker='o', s = 100, c='k')

    t, lum = tde()
    ax.scatter(
            t[0], lum[0], marker='o', edgecolor='k', facecolor='white',
            label="Rel. TDE", s=100)
    ax.text(t[0], lum[0]/3, 'SwiftJ1644+57', 
            verticalalignment='center', 
            horizontalalignment='center', fontsize=11)
    ax.annotate('', xy=(t[0]*1.5,lum[0]), xytext=(t[0],lum[0]),
            arrowprops=dict(
                facecolor='white', headwidth=10, width=1, headlength=7))
    ax.annotate('', xy=(t[0],lum[0]*3), xytext=(t[0],lum[0]),
            arrowprops=dict(
                facecolor='white', headwidth=10, width=1, headlength=7))


    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r"Time to Peak [days; rest frame]", fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlim(0.1,3E3)

    ax.set_ylim(1E36, 1E44)
    plt.tight_layout()
    plt.show()
    #plt.savefig("early_nu_lnu.png")
