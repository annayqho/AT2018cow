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
    inputf = direc + "/Soderberg2009Fig4.1"

    dat = Table.read(inputf, format='ascii')
    x = dat['col1']
    y = dat['col2']

    choose = np.logical_and(y < 1e50, x < 0.8)
    ax.scatter(
            x[choose], y[choose], marker='o', c='k', s=100,
            label="SNe Ibc")

    choose = np.logical_and(x > 0.5, x < 1.2)
    ax.scatter(
            x[choose], y[choose], marker='s', c='k', s=100,
            label="Rel. SNe")

    choose = y > 7e49
    ax.scatter(
            x[choose], y[choose], marker='o', edgecolors='k', s=100,
            facecolors='none', label=None)


    # AT2018cow
    ax.scatter(
            0.173, 2.2E49, 
            marker='*', s=300, facecolors='white', edgecolors='black')
    ax.text(
            0.173, 2.2E49/3, "AT2018cow", 
            fontsize=14, verticalalignment='bottom', horizontalalignment='center')


    # Swift J1644
    ax.scatter(
            0.5*1.2, 2.9E50 + 0.1*2.9E50, 
            marker='s', s=100, facecolors='white', edgecolors='black')

    
    # MAXI 140814A
    # ax.scatter(
    #          0.2, 3E49,
    #          marker='*', s=300, facecolors='white', edgecolors='black')
    # ax.text(
    #          0.2, 3E49*1.2, "VLASS 1210+49", 
    #          fontsize=14, verticalalignment='bottom', horizontalalignment='center')  


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
    direc = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
    inputf = direc + "/Soderberg2009Fig3.1"

    dat = Table.read(inputf, format='ascii')
    x = dat['col1']
    y = dat['col2']

    choose = np.logical_or(y < 9E27, x > 50)
    ax.scatter(
            x[choose], y[choose], marker='o', c='k', s=100,
            label="SNe Ibc")

    choose = np.logical_and(y > 8E27, x < 80)
    ax.scatter(
            x[choose], y[choose], marker='s', c='k', s=100,
            label="Rel. SNe")

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

    # MAXI 140814A
    # ax.scatter(
    #         4*365*6/5, 8e28,
    #         marker='*', s=300, facecolors='white', edgecolors='black')
    # ax.text(
    #         4*365*6/5, 9e28, "VLASS 1210+49", 
    #         fontsize=14, verticalalignment='bottom', horizontalalignment='center')

    # Swift J1644+57
    # d = Planck15.luminosity_distance(z=0.354).cgs.value
    # ax.scatter(
    #         22*80/5, 25*1e-3*1e-23*4*np.pi*d**2,
    #         marker='s', s=100, facecolors='white', edgecolors='black',
    #         label="Rel. TDE")

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
