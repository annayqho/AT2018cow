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
            facecolors='none', label="GRBs")


    # AT2018cow
    ax.scatter(
            0.173, 1.2E49, 
            marker='*', s=300, facecolors='white', edgecolors='black')
    ax.text(
            0.173, 2.2E49, "AT2018cow", 
            fontsize=14, verticalalignment='bottom', horizontalalignment='center')


    # Swift J1644
    ax.scatter(
            0.5*1.2, 2.9E50 + 0.1*2.9E50, 
            marker='s', s=100, facecolors='white', edgecolors='black',
            label="Rel. TDE")

    ax.legend(loc='upper left', fontsize=14)

    ax.set_xlim(0.03, 8)
    ax.set_ylim(1E45, 1E52)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlabel(
        "Blastwave Velocity $(\\Gamma \\beta)$",
        fontsize=14)
    ax.set_ylabel(
        "Energy (erg)", fontsize=14)

    

def peaklum(ax):
    direc = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
    inputf = direc + "/Soderberg2009Fig3.1"

    dat = Table.read(inputf, format='ascii')
    x = dat['col1']
    y = dat['col2']

    choose = np.logical_or(y < 9E27, x > 50)
    ax.scatter(
            x[choose], y[choose], marker='o', c='k', s=150,
            label="SNe Ibc")

    choose = np.logical_and(y > 8E27, x < 80)
    ax.scatter(
            x[choose], y[choose], marker='s', c='k', s=150,
            label="GRBs/Rel. SNe")

    # Lines
    rotangle = 55
    # v = 0.1c
    ax.plot([4.3,2300], [1E25,4.4E30], ls='--', c='k')
    ax.text(
            8, 2E25, "$R/\Delta t = 0.1c$", fontsize=14, rotation=rotangle,
            horizontalalignment='left', verticalalignment='bottom')

    # v = 0.01c
    ax.plot([43,6000], [1E25, 3E29], ls='--', c='k')
    ax.text(
            80, 2E25, "$R/\Delta t = 0.01c$", fontsize=14, rotation=rotangle,
            horizontalalignment='left', verticalalignment='bottom')

    # v = c
    ax.plot([1, 207], [7E25, 4e30], ls='--', c='k')
    ax.text(
            2, 1.8E26, "$R/\Delta t = c$", fontsize=14, rotation=rotangle,
            horizontalalignment='left', verticalalignment='bottom')

    # AT2018cow
    ax.scatter(
            22*100/5, 2E29, 
            marker='*', s=300, facecolors='white', edgecolors='black')
    ax.text(
            22*100/5, 3E29, "AT2018cow", 
            fontsize=14, verticalalignment='bottom', horizontalalignment='center')

    # MAXI 140814A
    # ax.scatter(
    #         4*365*6/5, 8e28,
    #         marker='*', s=300, facecolors='white', edgecolors='black')
    # ax.text(
    #         4*365*6/5, 9e28, "MAXI 140814A", 
    #         fontsize=14, verticalalignment='bottom', horizontalalignment='center')

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
vele(axarr[0])
peaklum(axarr[1])
plt.tight_layout()

#plt.show()
plt.savefig("chevalier_diagram.png")
