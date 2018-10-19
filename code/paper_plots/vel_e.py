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

squaresize = 50

def vele(ax):
    # only use this to plot the GRBs
    direc = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
    inputf = direc + "/Soderberg2009Fig4.1"

    dat = Table.read(inputf, format='ascii')
    x = dat['col1']
    y = dat['col2']

    choose = np.logical_and(x > 1, y > 3E49)
    ax.scatter(
            x[choose], y[choose], marker='x', c='k', s=50, label="GRBs")

    # Swift TDE
    # just use one epoch
    # epoch 3: day 18
    v = 1.7
    E = 2.7E51
    ax.scatter(
            v, E, marker='o', edgecolor='k', facecolor='k', s=50, label="TDE")
    ax.text(
            v/1.1, E, "Swift J1644", fontsize=12,
            horizontalalignment='right',
            verticalalignment='top')

    # SN 2003L
    v = 0.12
    E = 6.9E47
    ax.scatter(
            v, E, marker='+', c='k', s=100, label="SNe Ibc")
    ax.text(
            v*1.1, E, "2003L", fontsize=12,
            horizontalalignment='left',
            verticalalignment='center')

    # SN 2007bg
    v = 0.19
    E = 2.5E48
    ax.scatter(
            v, E, marker='+', c='k', s=100, label=None)
    ax.text(
            v, E*1.3, "2007bg", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='center')

    # SN 2003bg
    v = 0.11
    E = 8.7E47
    ax.scatter(
            v, E, marker='+', c='k', s=100, label=None)
    ax.text(
            v, E*1.3, "2003bg", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='center')

    # SN 1998bw
    v = 1.26
    E = 4.8E48
    ax.scatter(
            v, E, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label="LLGRB-SNe")
    ax.text(
            v, E*1.3, "1998bw", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='center')

    # SN 2009bb
    v = 0.71
    E = 2.9E48
    ax.scatter(
            v, E, marker='+', c='k', s=100)
    ax.text(
            v, E*1.3, "2009bb", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='center')

    # SN 2006aj
    v = 2.1
    E = 7.4E47
    ax.scatter(
            v, E, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label=None)
    ax.text(
            v, E*1.3, "2006aj", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='center')

    # SN 2010bh
    v = 0.33
    E = 9.0E47
    ax.scatter(
            v, E, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label=None)
    ax.text(
            v, E*1.3, "2010bh", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='center')


    # SN 1988Z
    v = 0.01
    E = 2.0E48
    ax.scatter(
            v, E, marker='o', edgecolor='k', s=100,
            facecolor='white', label="SNe II")
    ax.text(
            v, E*1.3, "88Z", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='center')

    # SN 1979C
    v = 0.016
    E = 1E48
    ax.scatter(
            v, E, marker='o', edgecolor='k', s=100,
            facecolor='white', label=None)
    ax.text(
            v, E*1.3, "79C", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='center')




    # AT2018cow
    ax.scatter(
            0.173, 2.2E49, 
            marker='*', s=300, facecolors='black', edgecolors='black')
    ax.text(
            0.173, 2.2E49*1.3, "AT2018cow", 
            fontsize=14, 
            verticalalignment='bottom', horizontalalignment='center')


    ax.set_xscale('log')
    ax.set_yscale('log')

    # make a twin axis
    ax2 = ax.twinx()
    ax2.set_ylabel(
            "Energy (erg) $= U_B/\epsilon_B$, $\qquad \epsilon_B=0.01$", 
            fontsize=14, rotation=270, labelpad=15.0)
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
            "Energy (erg) $= U_B/\epsilon_B$, \qquad $\epsilon_B=0.33$", 
            fontsize=14)
    ax.set_xlim(0.007, 8)
    ax.set_ylim(4E47, 4E51)

    ax.legend(loc='upper left', fontsize=10)

fig,ax = plt.subplots(1,1, figsize=(5,5))
vele(ax)
plt.tight_layout()

#plt.show()
plt.savefig("vel_e.png")
