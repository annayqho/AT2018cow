""" Plot nu_peak over time for various transients, including AT2018cow """


import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)


def sn2009bb(ax):
    dt = [20, 52, 81, 145]
    nu_peak = [6, 3, 1, 0.8]
    ax.scatter(dt, nu_peak, c='grey', marker='o')
    ax.plot(dt, nu_peak, c='grey', lw=2.0)
    ax.text(
            dt[-1], nu_peak[-1], "2009bb", 
            horizontalalignment='left', fontsize=12,
            verticalalignment='center')



def extra():
    # I don't know the peak frequency

    plt.errorbar(5, 3.15, yerr=1.75, fmt='.', c='grey') 
    plt.text(5, 3.15, "1998bw", fontsize=12)


    # Some IIn's and stuff
    # http://adsabs.harvard.edu/abs/2015aska.confE..60P

    # SN 1993J
    # http://iopscience.iop.org/article/10.1086/523258/pdf

    dt = [16.4, 20.4, 25.4, 35.3, 45.3, 58.2, 92.9]
    nu_peak = [10, 8.13, 6.18, 4.72, 3.91, 3.184, 2.235]
    plt.scatter(dt, nu_peak, c='grey', marker='o')
    plt.plot(dt, nu_peak, c='grey', marker='o')
    plt.text(dt[-1], nu_peak[-1], "SN 2011dh (IIb)", fontsize=12)


    # https://arxiv.org/pdf/astro-ph/0512413.pdf
    # Table 33

    dt = 365*np.array([2.6, 1.2, 0.4, 3.8, 0.0093, 2.5, 0.030, 0.20, 0.27, 0.068])
    nu_peak = np.array([4.9, 5.0, 5.0, 4.9, 0.84, 4.9, 4.9, 5.0, 4.9, 8.5])
    names = np.array(
            ['1978K', '1979C', '1980K', '1986J', '1987A', '1988Z', '1998bw',
             '2001ig', '2003bg', '2004cc'])
    plt.scatter(dt, nu_peak, c='grey', marker='s')
    for ii,val in enumerate(dt):
        plt.text(val, nu_peak[ii], names[ii], fontsize=10)

    # 2003bg
    plt.scatter([23, 351], [25, 2.2], c='grey', marker='s')
    plt.plot([23, 351], [25, 2.2], c='grey', marker='s')
    plt.text(351, 2.2, "2003bg", fontsize=14)


    # 2003L
    # http://iopscience.iop.org/article/10.1086/427649/pdf
    plt.scatter([30], [22.5], c='grey', marker='s')
    plt.text(30, 22.5, "2003L")




fig, ax = plt.subplots(1, 1, figsize=(10,10))

# AT 2018cow
ax.scatter(22, 92, c='black', marker='s')
ax.text(
        22, 92, "AT2018cow", fontsize=14,
        verticalalignment='center')

# plot

ax.set_ylabel("$\\nu_p$ [GHz]", fontsize=16)
ax.set_xlabel("$\\Delta t$ [days]", fontsize=16)
ax.tick_params(axis='both', labelsize=14)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(2, 2000)
ax.set_ylim(0.5, 200)
plt.tight_layout()

plt.show()
#plt.savefig("nupeak_evolution.png")
