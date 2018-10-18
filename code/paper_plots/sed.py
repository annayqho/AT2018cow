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
from get_optical import get_opt_sed
from xray_lc import get_xrt, get_nustar
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
    t_x, f_x, ef_x = get_xrt()
    lum_x = f_x * 4 * np.pi * d**2
    freq.append(1E17)
    nufnu.append(np.interp(day, t_x, lum_x))
    band.append('x')

    # radio
    nu_rad, fnu_rad = get_spectrum(day)
    for ii,val in enumerate(nu_rad):
        freq.append(val*1E9)
        nufnu.append(val*1E9*mjy_to_lum(fnu_rad[ii]))
        band.append('r')

    # optical
    nu_opt, fnu_opt = get_opt_sed(day)
    for ii,val in enumerate(nu_opt):
        freq.append(val)
        nufnu.append(val*mjy_to_lum(fnu_opt[ii]))
        band.append('o')

    freq = np.array(freq)
    nufnu = np.array(nufnu)
    band = np.array(band)
    order = np.argsort(freq)
    freq = freq[order]
    nufnu = nufnu[order]
    band = band[order]

    return band, freq, nufnu


def ext_ind(b, x, y, day, col):
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
        yvals = yr[-1] * (xvals/xr[-1])**(plot_ind)
        plt.plot(xvals, yvals, c=col, ls='--')


def plot_ind(b, x, y, day, col):
    """ Extrapolate spectral index for a certain day """
    days = np.array([22])
    ind = np.array([-0.7])
    new_ind = 1 + ind
    if sum(days==day) > 0:
        plot_ind = new_ind[days==day]
        # plot from the last radio point
        xr = x[b == 'r']
        yr = y[b == 'r']
        xvals = np.linspace(xr[-1], 2E16)
        yvals = yr[-1] * (xvals/xr[-1])**(plot_ind)
        plt.plot(xvals, yvals, c='k', ls=':', lw=1)


def plot_xray_spindex(x, y, col):
    """ Dashed line indicating X-ray spectral index """
    Gamma = 1.54
    # f_nu \propto nu^{-alpha}
    # Gamma = photon index = alpha + 1
    # 0.3 -- 10 keV
    # 7.2E16 Hz -- 2.4E18 Hz
    xvals = np.linspace(2.4E16, 2.4E18)
    yvals = y * (xvals/x)**(0.46)
    plt.plot(xvals, yvals, ls='-', c=col)


def plot_xray(flux, c):
    """ OK so the approach is to take the integrated flux across 0.3-10 keV,
    and use that and the geometric mean of the frequency and the spectral index
    to solve for the normalization coefficient for your spectrum

    Parameters
    ----------
    x: frequency
    y: flux
    c: color of the point/line

    2.4E18 Hz corresponds to 10 keV
    7.2E16 Hz corresponds to 0.3 keV
    """
    nu0 = np.sqrt(0.3*10) * (2.4E18/10) 
    alpha = 0.54
    nu2 = 2.4E18
    nu1 = 7.2E16
    A = flux*(1-alpha) / (nu0**alpha * (nu2**(1-alpha) - nu1**(1-alpha)))
    xplot = np.linspace(7.2E16, 2.4E18, 1000)
    yplot = A*(xplot/nu0)**(-alpha)
    plt.plot(xplot, xplot*yplot, c=c, ls='-')
    if c=='white':
        plt.plot(xplot, xplot*yplot, c='k', ls='--')




# set up the plot
fig = plt.figure(figsize=(6,4))

# the days that have ALMA data:
# 14, 22, 24

days = [5, 14, 22]
#cols = ['#1b9e77', '#d95f02', '#7570b3']
#cols = ['#440154', '#5ec962', '#fde725']
cols = ['#f98e09', '#bc3754', '#57106e']
#cols = ['k', 'white', 'lightgrey']
markers = ['o', 's', 'D']

for ii,day in enumerate(days):
    b, x, y = get_sed(day)
    choose = b != 'x'
    plt.scatter(
            x[choose], y[choose], edgecolor='k', facecolor=cols[ii], 
            label=r"Day %s" %(day), marker=markers[ii], lw=0.5)
    # For Day 14, choose the highest-frequency radio point and draw
    # a line connecting 340E9 Hz to 1.5E14 Hz
    xplt = np.logspace(np.log10(340E9), np.log10(1.5E14))
    yplt = 4E40*(xplt/230E9)**(1-0.75)
    plt.plot(xplt, yplt, ls='--', c='k', lw=0.5)
    ind = int(len(xplt)/2)
    plt.text(xplt[ind], yplt[ind], '$\propto \\nu^{-0.75}$', fontsize=11,
            horizontalalignment='center', verticalalignment='bottom')
    yplt = 4E40*(xplt/230E9)**(1-1)
    plt.plot(xplt, yplt, ls='--', c='k', lw=0.5)
    plt.text(xplt[ind], yplt[ind]/1.5, '$\propto \\nu^{-1}$', fontsize=11,
            horizontalalignment='center', verticalalignment='top')
    choose = b == 'x'
    plot_xray(y[choose][0], cols[ii])
    if ii == 2:
        plt.text(
                x[choose][0], y[choose][0], 
                '$\propto \\nu^{0.54}$', fontsize=12)



# Band 9 point
nu = 671E9
f = 31.5
lum = mjy_to_lum(f)
plt.errorbar(
        nu, nu*lum, yerr=0.2*lum,
        fmt='*', mec='k', mfc='white', ms=13, mew=0.5)
        #label=r"Day 24")

nu_opt, fnu_opt = get_opt_sed(24)
lum_opt = nu_opt*mjy_to_lum(fnu_opt)
plt.errorbar(
        nu_opt, lum_opt, 
        fmt='*', mec='k', mfc='white', label=r"Day 24", ms=13, mew=0.5)

# Plot the NuSTAR data
dt, f1, ef1, f2, ef2, f3, ef3, f4, ef4 = get_nustar()
cols = ['#08519c', '#3182bd', '#6baed6', '#bdd7e7']

def plot_nustar(f, ef, minf, maxf):
    center = (minf+maxf)/2
    ec = maxf-center
    lum = f * 10**(-12) * 4 * np.pi * d**2
    elum = ef * 10**(-12) * 4 * np.pi * d**2
    for ii in np.arange(len(lum)):
        plt.errorbar(
                center, lum[ii], 
                xerr=ec, yerr=elum[ii], c=cols[ii])

#plot_nustar(f1, ef1, 7.25e17, 2.42e18) # 3-10 keV
#plot_nustar(f2, ef2, 2.41798926e18, 4.84E18) # 10-20 keV
#plot_nustar(f3, ef3, 4.84E18, 9.67E18) # 20-40 keV
#plot_nustar(f4, ef4, 9.67E18, 1.93E19) # 40-80 keV

# Model for Epoch 2 is a power-law:
norm = 9E-4 # photons/keV/cm2/s at 1 keV
f = norm * 10**(-12)  # maybe in units 1E-12 erg/cm2/s
alpha = 1.39
nu = np.linspace(7.25E17, 1.93E19, 1000)
energy = nu * 6.626E-27
y = f * energy**(-alpha) * 4 * np.pi * d**2
#plt.plot(nu, y, ls=':', c='k')

plt.xscale('log')
plt.yscale('log')
plt.xlabel("Frequency [Hz]", fontsize=16)
plt.ylabel(r"$\nu L_\nu$ (erg/s)", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12, loc='upper left')

plt.tight_layout()

plt.savefig("sed.png")
#plt.show()
