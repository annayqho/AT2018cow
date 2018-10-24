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


def vel_lines(ax, x, v):
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
    ax.plot(xvals, yvals, ls='--', c='k', lw=0.5)
    rotangle = 70
    ax.text(
            x, 7E27, "$R/\Delta t = %sc$" %v, 
            fontsize=10, rotation=rotangle,
            horizontalalignment='center', verticalalignment='top')
    return yvals


def mdot_curves(ax, x, mdotv):
    # Mdot divided by v
    xvals = np.linspace(1,3000)
    logy = (19/4) * np.log10(0.002) - (19/4)*np.log10(mdotv) + \
            (2*19/4)*np.log10(xvals) 
    yvals = 10**logy
    ax.plot(xvals, yvals, ls=':', c='k', lw=0.5)
    rotangle = 84
    ax.text(
            x, 5E29, 
            "$\dot{M}/v = 10^{%s}$" %int(np.log10(mdotv)), 
            fontsize=10, rotation=rotangle,
            horizontalalignment='left', verticalalignment='top')
    return yvals


def density_curves(ax, x, ne):
    """ 
    Parameters
    ----------
    ax: axis to plot the stuff on
    x: the value of (dt / 1 d) * (nu_p / 5 GHz)
    v: density in units of parts per cm cubed
    """
    xvals = np.linspace(1,3000)
    logy = (19/22) * np.log10(79) + 26 - (19/22)*np.log10(ne) + \
            (19/22) * 2 * np.log10(xvals**2 / 22) # divide out 22 days
    yvals = 10**logy
    ax.plot(xvals, yvals, ls=':', c='k', lw=0.5)
    rotangle = 75 
    ax.text(
            x, 5E29, "$n_e = 10^{%s} \mathrm{cm}^{-3}$" %int(np.log10(ne)), 
            fontsize=10, rotation=rotangle,
            horizontalalignment='left', verticalalignment='top')
    return yvals


def lumtnu(ax):
    # ASASSN14li
    tnu = (143)*(8.20/5)
    lpeak = 1.8E28
    ax.scatter(
            tnu, lpeak, marker='o', edgecolor='k', facecolor='black', s=100,
            label='TDE')
    ax.text(
            tnu, lpeak/1.3, "ASASSN14li", fontsize=12,
            verticalalignment='top',
            horizontalalignment='center')

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
            label=None)
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
    y = vel_lines(ax, 7, 1)
    y = vel_lines(ax, 70, 0.1)
    y = vel_lines(ax, 700, 0.01)

    # AT2018cow
    ax.scatter(
            22*100/5, 4.4E29, marker='*', s=300, 
            facecolors='black', edgecolors='black')
    ax.text(
            22*100/5, 5E29, "AT2018cow", fontsize=14, 
            verticalalignment='bottom', 
            horizontalalignment='center')


    ax.set_xlim(2, 3000)
    ax.set_ylim(5E26, 1E30)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xlabel(
        "$(\Delta t/1\,\mathrm{day})(\\nu_p/5\,\mathrm{GHz})$",
        fontsize=14)

    
fig,axarr = plt.subplots(1,2, figsize=(8,5))
ax = axarr[0]
lumtnu(ax)
y = density_curves(ax, 7, 1E-1)
y = density_curves(ax, 70, 1E3)
y = density_curves(ax, 700, 1E7)
ax.set_ylabel("Peak Radio Luminosity ($\mathrm{erg\,s^{-1}\,Hz^{-1}}$)",
    fontsize=14)
ax = axarr[1]
lumtnu(ax)
y = mdot_curves(ax, 18, 1E-8)
y = mdot_curves(ax, 180, 1E-6)
y = mdot_curves(ax, 1800, 1E-4)
y = mdot_curves(ax, 1.8, 1E-2)
ax.get_yaxis().set_visible(False)
ax.legend(loc='lower center', ncol=2)#, columnspacing=0.1)
#y = mdot_curves(ax, 700, 1E1)

plt.tight_layout()


#plt.show()
plt.savefig("lum_tnu.png")
