""" Peak lum phase space from Soderberg (2009) with AT2018cow
and Swift J1644 on top

Combine with the Chevalier diagram

These are two classic figures
"""

from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.cosmology import Planck15


def chev_lines(ax, x, v):
    """ Equation 16 from Chevalier 1998 
    
    Parameters
    ----------
    ax: axis to plot the stuff on
    x: the value of (dt / 1 d) * (nu_p / 5 GHz)
    v: velocity in units of c
    """
    logy = (36/17)*\
            ((17*26/36) + \
            np.log10(v*3e10/(10*3.1E9)) + \
            np.log10(x))
    y = 10**logy
    ax.plot(x, y, ls='--', c='k')
    rotangle = 60
    ax.text(
            x[5300], y[5500], "$R/\Delta t = %sc$" %v, 
            fontsize=14, rotation=rotangle,
            horizontalalignment='center', verticalalignment='top')
    return y


def vele(ax):
    direc = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"

    # SN 2003L
    v = 0.16
    E = 6.9E47
    ax.scatter(
            v, E, marker='o', edgecolor='k', s=100,
            facecolor='none', label="SNe Ibc")
    ax.text(
            v*1.05, E, "2003L", fontsize=12)

    # SN 2007bg
    v = 0.49
    E = 8.2E47
    ax.scatter(
            v, E, marker='o', edgecolor='k', s=100,
            facecolor='none', label="SNe Ibc")
    ax.text(
            v*1.05, E, "2007bg", fontsize=12)

    #ax.scatter(
    #        x[choose], y[choose], marker='s', c='k', s=100,
    #        label="Rel. SNe")

    #choose = y > 7e49
    #ax.scatter(
    #        x[choose], y[choose], marker='o', edgecolors='k', s=100,
    #        facecolors='none', label=None)


    # Mark 2003bg
    #ax.text(
    #    0.24, 1E48, "SN2003bg", fontsize=12)

    # AT2018cow
    ax.scatter(
            0.173, 2.2E49, 
            marker='*', s=300, facecolors='white', edgecolors='black')
    ax.text(
            0.173, 2.2E49/3, "AT2018cow", 
            fontsize=14, verticalalignment='bottom', horizontalalignment='center')


    ax.set_xscale('log')
    ax.set_yscale('log')

    # make a twin axis
    ax2 = ax.twinx()
    ax2.set_ylabel(
            "Energy (erg): $\epsilon_B=0.01, \epsilon_e=0.1$", fontsize=14,
            rotation=270, labelpad=15.0)
    y_f = lambda y_i: y_i*33
    ymin, ymax = ax.get_ylim()
    ax2.set_ylim((y_f(ymin), y_f(ymax)))
    ax2.plot([],[])
    ax2.set_yscale('log')

    ax.tick_params(axis='both', labelsize=14)
    ax2.tick_params(axis='both', labelsize=14)
    ax.set_xlabel(
        "Blastwave Velocity $(\\Gamma \\beta)$",
        fontsize=14)
    ax.set_ylabel(
            "Energy (erg): $\epsilon_B=\epsilon_e=0.33$", fontsize=14)
    ax.set_xlim(0.03, 2)
    ax.set_ylim(1E45, 8E49)


def peaklum(ax):
    # 2003L
    tnu = (30)*(22.5/5)
    lpeak = 3.3E28
    ax.scatter(
            tnu, lpeak, marker='o', edgecolor='k', s=100,
            facecolor='none',
            label="SNe Ibc")
    ax.text(
            tnu*1.05, lpeak, "2003L", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='left')

    # 2007bg
    tnu = (55.9)*(8.46/5)
    lpeak = 4.1E28
    ax.scatter(
            tnu, lpeak, marker='o', edgecolor='k', s=100,
            facecolor='none',
            label="SNe Ibc")
    ax.text(
            tnu*1.05, lpeak, "2007bg", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='right')

    #ax.scatter(
    #        x[choose], y[choose], marker='s', c='k', s=100,
    #        label="Rel. SNe")

    # Lines
    x = np.logspace(0, 4, 10000)
    y = chev_lines(ax, x, 0.1)
    y = chev_lines(ax, x, 1)
    y = chev_lines(ax, x, 0.01)

    # AT2018cow
    ax.scatter(
            22*100/5, 2E29, marker='*', s=300, 
            facecolors='white', edgecolors='black')
    ax.text(
            22*100/5, 3E29, "AT2018cow", fontsize=14, 
            verticalalignment='bottom', 
            horizontalalignment='center')


    # Mark 2003bg
    # ax.text(
    #         23*25/5, 5.5E28, "2003bg", fontsize=12,
    #         verticalalignment='top',
    #         horizontalalignment='left')

    # Mark SN2007bg
    # luminosity reaches 1E29 567 days after the explosion, at 8.46 GHz
    # ax.scatter(
    #         55.9*8.46/5, 4.1E28, marker='o', edgecolor='k', s=100,
    #         facecolor='none')
    # ax.text(
    #         55.9*8.46/5, 1.2*4.1E28, "SN2007bg", fontsize=12,
    #         verticalalignment='bottom',
    #         horizontalalignment='center')

    ax.legend(loc='upper left', fontsize=14)

    ax.set_xlim(1, 6000)
    ax.set_ylim(1E25, 5E30)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlabel(
        "$(\Delta t/1\,\mathrm{day})(\\nu_p/5\,\mathrm{GHz})$",
        fontsize=14)
    ax.set_ylabel("Peak Radio Luminosity ($\mathrm{erg\,s^{-1}\,Hz^{-1}}$)",
        fontsize=14)

    
fig,axarr = plt.subplots(1,2, figsize=(10,5))
vele(axarr[1])
peaklum(axarr[0])
plt.tight_layout(w_pad=4.0)

#plt.show()
plt.savefig("chevalier_diagram.png")
