""" 
Plot the spectrum from radio to X-ray frequencies 

Choose Days 10 (June 26), 13/14 (June 29/30), and 22 (July 8),
because these have the best sampling.
"""

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import numpy as np
from astropy.cosmology import Planck15
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/code")
from get_optical import get_opt_data
from xray_lc import get_xray
from get_radio import get_spectrum, get_data_all

d = Planck15.luminosity_distance(z=0.014).cgs.value


def get_date(dt):
    """ Convert dt to date UT """
    day = dt + 16
    month = "June"
    if day > 30:
        month = "July"
        day = day - 30
    return "%s %s" %(day, month)


def mjy_to_lum(f):
    """ Convert flux from mJy to erg/s/Hz """
    lum = f * 1e-3 * 1e-23 * 4 * np.pi * d**2
    return lum


def get_sed(day):
    """ Return freq and nufnu for an individual day """
    freq = []
    nufnu = []
    band = []

    # X rays
    t_x, f_x, ef_x = get_xray()
    lum_x = f_x * 4 * np.pi * d**2
    freq.append(1E17)
    nufnu.append(np.interp(day, t_x, lum_x))
    band.append('x')
    
    # radio
    nu_rad, fnu_rad = get_spectrum(day)
    for ii,val in enumerate(nu_rad):
        freq.append(val*1E9)
        nufnu.append(val * 1E9 * mjy_to_lum(fnu_rad[ii]))
        band.append('r')

    # optical
    # goes from 1,000 -- 20,000 A
    # which is 10^14 -- 10^15 Hz
    # I need to ask Dan for the SED
    t_o, lum_o = get_opt_data()
    freq.append(1E14)
    nufnu.append(np.interp(day, t_o, lum_o))
    band.append('o')

    freq = np.array(freq)
    nufnu = np.array(nufnu)
    band = np.array(band)
    order = np.argsort(freq)
    freq = freq[order]
    nufnu = nufnu[order]
    band = band[order]

    return band, freq, nufnu


def plot_ind(b, x, y, day, col):
    """ Extrapolate spectral index for a certain day """
    days = np.array([10, 14, 22])
    ind = np.array([-1.8, -0.8, -1.2])
    new_ind = 1 + ind
    if sum(days==day) > 0:
        plot_ind = new_ind[days==day]
        # plot from the last radio point
        xr = x[b == 'r']
        yr = y[b == 'r']
        xvals = np.linspace(xr[-1], 1E18)
        print(plot_ind)
        yvals = yr[-1] * (xvals/xr[-1])**(plot_ind)
        plt.plot(xvals, yvals, c=col, ls='--')


def plot_xray(x, y, col):
    """ Dashed line indicating X-ray spectral index """
    Gamma = 1.54
    # f_nu \propto nu^{-alpha}
    # Gamma = photon index = alpha + 1
    # 0.1 -- 10 keV
    # 2.4E16 Hz -- 2.4E18 Hz
    xvals = np.linspace(2.4E16, 2.4E18)
    yvals = y[-1] * (xvals/x[-1])**(0.46)
    plt.plot(xvals, yvals, ls='-', c=col)


# set up the plot
fig = plt.figure(figsize=(6,4))


# the days that have ALMA data:
# 14, 22, 24

days = [5, 14, 22]
cols = ['k', 'grey', 'lightgrey']
markers = ['o', 's', 'D']

for ii,day in enumerate(days):
    b, x, y = get_sed(day)
    plt.scatter(
            x, y, edgecolor='k', facecolor=cols[ii], 
            label=r"Day %s" %(day), marker=markers[ii], lw=0.5)
    #choose = b == 'r'
    #plt.plot(x[choose], y[choose], c=cols[ii])
    plot_xray(x, y, cols[ii])
    plot_ind(b, x, y, day, cols[ii])

# Band 9 point
nu = 671E9
f = 31.5
lum = mjy_to_lum(f)
plt.errorbar(
        nu, nu*lum, yerr=0.2*lum,
        fmt='*', mec='k', mfc='white', label=r"Day 24", ms=13, mew=0.5)

plt.xscale('log')
plt.yscale('log')
plt.xlabel("Frequency [Hz]", fontsize=16)
plt.ylabel(r"$\nu L_\nu$ (erg/s)", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()

plt.tight_layout()

plt.savefig("sed.png")
#plt.show()
