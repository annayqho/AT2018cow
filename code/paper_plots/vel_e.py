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


def chev_lines(ax, x, v):
    """ Equation 16 from Chevalier 1998 
    
    Parameters
    ----------
    ax: axis to plot the stuff on
    x: the value of (dt / 1 d) * (nu_p / 5 GHz)
    v: velocity in units of c
    """
    xvals = np.linspace(1,3000)
    logy = (26) + \
            (19/9) * np.log10(v) + \
            (19/9) * np.log10(xvals)
    yvals = 10**logy
    ax.plot(xvals, yvals, ls='--', c='k')
    rotangle = 70
    ax.text(
            x, 7E27, "$R/\Delta t = %sc$" %v, 
            fontsize=14, rotation=rotangle,
            horizontalalignment='center', verticalalignment='top')
    return yvals


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

    ax.legend(loc='upper left', fontsize=14)



def peaklum(ax):
    # 88Z
    tnu = (1253)*(5/5)
    lpeak = 2.2E28
    ax.scatter(
            tnu, lpeak, marker='o', edgecolor='k', facecolor='white', s=100,
            label='SN II')
    ax.text(
            tnu, lpeak/1.2, "88Z", fontsize=12,
            verticalalignment='top',
            horizontalalignment='right')

    # 79C 
    tnu = (1400)*(1.4/5)
    lpeak = 4.3E27
    ax.scatter(
            tnu, lpeak, marker='o', edgecolor='k', facecolor='white', s=100,
            label='SN II')
    ax.text(
            tnu, lpeak/1.2, "79C", fontsize=12,
            verticalalignment='top',
            horizontalalignment='right')

    # 2003L
    tnu = (30)*(22.5/5)
    lpeak = 3.3E28
    ax.scatter(
            tnu, lpeak, marker='+', c='k', s=100,
            label="SNe Ibc")
    ax.text(
            tnu, lpeak/1.2, "2003L", fontsize=12,
            verticalalignment='top',
            horizontalalignment='right')

    # 11qcj
    tnu = (10)*(5/5)
    lpeak = 7E28
    ax.scatter(
            tnu, lpeak, marker='+', c='k', s=100,
            label=None)
    ax.text(
            tnu/1.2, lpeak, "11qcj", fontsize=12,
            verticalalignment='center',
            horizontalalignment='right')


    # 2007bg
    tnu = (55.9)*(8.46/5)
    lpeak = 4.1E28
    ax.scatter(
            tnu, lpeak, marker='+', c='k', s=100, label=None)
    ax.text(
            tnu, lpeak*1.2, "2007bg", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='center')

    # SN 2003bg
    tnu = (35)*(22.5/5)
    lpeak = 3.9E28
    ax.scatter(
            tnu, lpeak, marker='+', c='k', s=100, label=None)
    ax.text(
            tnu*1.2, lpeak, "2003bg", fontsize=12,
            verticalalignment='center',
            horizontalalignment='left')

    # SN 1998bw
    tnu = (10)*(10/5)
    lpeak = 8.2E28
    ax.scatter(
            tnu, lpeak, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label="LLGRB-SNe")
    ax.text(
            tnu, lpeak*1.2, "1998bw", fontsize=12,
            verticalalignment='bottom',
            horizontalalignment='center')

    # SN 2009bb
    tnu = (20)*(6/5)
    lpeak = 3.6E28
    ax.scatter(
            tnu, lpeak, marker='+', c='k', s=100, label=None)
    ax.text(
            tnu/1.2, lpeak, "2009bb", fontsize=12,
            verticalalignment='center',
            horizontalalignment='right')

    # SN 2006aj
    tnu = (5)*(4/5)
    lpeak = 8.3E27
    ax.scatter(
            tnu, lpeak, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label=None)
    ax.text(
            tnu, lpeak/1.5, "2006aj", fontsize=12,
            verticalalignment='top',
            horizontalalignment='center')

    # SN 2010bh
    tnu = (30)*(5/5)
    lpeak = 1.2E28
    ax.scatter(
            tnu, lpeak, marker='s', edgecolor='k', s=squaresize,
            facecolor='k', label=None)
    ax.text(
            tnu, lpeak/1.5, "2010bh", fontsize=12,
            verticalalignment='top',
            horizontalalignment='center')

    # Lines
    y = chev_lines(ax, 7, 1)
    y = chev_lines(ax, 70, 0.1)
    y = chev_lines(ax, 700, 0.01)

    # AT2018cow
    ax.scatter(
            22*100/5, 2E29, marker='*', s=300, 
            facecolors='black', edgecolors='black')
    ax.text(
            22*100/5, 3E29, "AT2018cow", fontsize=14, 
            verticalalignment='bottom', 
            horizontalalignment='center')


    ax.set_xlim(2, 3000)
    ax.set_ylim(1E27, 1E30)
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

plt.show()
#plt.savefig("chevalier_diagram.png")
