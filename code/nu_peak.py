""" Plot nu_peak over time for various transients, including AT2018cow """


import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.table import Table
from plot_spec import fit_self_abs
from get_1993J import format_1993J


def get_peak(ax, freq_raw, flux_raw, flux_raw_err=None):
    """ 
    Estimate the peak frequency taking two power laws.

    If there's only one point, you can't do anything.
    """
    if len(freq_raw) < 2:
        print("Error: need at least 2 points")

    else:
        # sort by increasing frequency
        order = np.argsort(freq_raw)
        freq = freq_raw[order]
        flux = flux_raw[order]

        # this is the grid you will interpolate onto
        xgrid = np.linspace(min(freq), max(freq), 1000)

        # the peak flux is the first point of the second half
        # and the last point of the first half
        ind = np.argmax(flux)

        # if the peak is the first or last point, then it's a limit
        if ind == 0:
            # if the peak is first
            return (freq[0], True)
        elif ind == len(freq)-1:
            # if the peak is last
            return (freq[-1], True)
        else:
            # then you can actually do a fit
            m,b = np.polyfit(
                    np.log10(freq[0:2]), np.log10(flux[0:2]), deg=1)
            yfit_up = 10**(m*np.log10(xgrid)+b)
            #ax.plot(xgrid, yfit_up, ls=':', c='k')

            # plot and fit a nu*something line
            m,b = np.polyfit(np.log10(freq[2:]), np.log10(flux[2:]), deg=1)
            yfit_down = 10**(m*np.log10(xgrid)+b)
            #ax.plot(xgrid, yfit_down, ls=':', c='k')

            nupeak = xgrid[np.argmin(np.abs(yfit_up-yfit_down))]
            #ax.axvline(x=nupeak, ls='--')
            return (nupeak, False)


def sn2009bb(ax, color, label):
    """ Plot for SN 2009bb """
    dt = [20, 52, 81, 145]
    nu_peak = [6, 3, 1, 0.8]
    ax.scatter(dt[0:2], nu_peak[0:2], c=color, marker='o')
    ax.plot(dt[0:2], nu_peak[0:2], c=color, lw=1.0, label=label)
    ax.text(
            dt[0], nu_peak[0], "SN2009bb", 
            horizontalalignment='right', fontsize=14,
            verticalalignment='center')


def at2018cow(ax):
    """ Plot for AT 2018cow """

    dt = np.array([5, 7, 8, 10, 13.5, 22])
    islim = np.array([1, 0, 0, 0, 0, 0])
    nu = np.array([231.5, 281.15, 223.5, 144, 111, 92])

    choose = islim == 1
    ax.scatter(dt[choose], nu[choose], c='black', marker='^')

    choose = islim == 0
    ax.scatter(dt[choose], nu[choose], c='black', marker='s')
    ax.plot(dt[choose], nu[choose], c='black', lw=2.0)
    ax.text(
            dt[choose][-1], nu[choose][-1], "AT2018cow", fontsize=14,
            verticalalignment='center')


def sn1993J(ax, color):
    """ Plot for SN 1993J """
    dt, freq, flux = format_1993J()
    order = np.argsort(dt)
    choose = flux[order] > 0
    dt = dt[order][choose]
    freq = freq[order][choose]
    flux = flux[order][choose]
    ax.scatter(11, 87, marker='v', c=color) 
    ax.scatter(32, 25, marker='v', c=color)
    ax.scatter(75, 15, marker='v', c=color)
    ax.scatter(95, 15, marker='v', c=color)
    ax.plot([11, 32, 75, 95], [87, 25, 15, 15], ls='--', c=color)
    ax.text(
            11, 87, "SN1993J", fontsize=14, horizontalalignment='right')


def sn1998bw(ax, color, label):
    """ Plot for SN 1998bw 
    http://adsabs.harvard.edu/abs/1999A%26AS..138..467W
    """
    data_dir = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
    dat = Table.read(
            "%s/radio_compilations/1998bw.dat" %data_dir, 
            format='ascii', delimiter='&')
    dt = np.array(dat['dt'])
    flux = np.array(dat['flux'])
    freq = np.array(dat['freq'])

    dt_keep = []
    nu = []
    islim = []

    for day in dt:
        choose = dt == day
        if sum(choose) > 1:
            dt_keep.append(day)
            # fig, ax = plt.subplots(1,1)
            nupeak_val,islim_val = get_peak(
                    ax, freq[choose], flux[choose])
            # ax.scatter(
            #         freq[choose][np.argsort(freq[choose])], 
            #         flux[choose][np.argsort(freq[choose])])
            # ax.set_xscale('log')
            # ax.set_yscale('log')
            # plt.savefig("%s.png" %day)
            # plt.close()
            nu.append(nupeak_val)
            islim.append(islim_val)

    dt_keep = np.array(dt_keep)
    nu = np.array(nu)
    islim = np.array(islim)

    #choose = islim == True
    #ax.scatter(dt_keep[choose], nu[choose], c='grey', marker='^')

    choose = islim == False
    order = np.argsort(dt_keep[choose])
    ax.scatter(dt_keep[choose][order], nu[choose][order], c=color, marker='s')
    ax.plot(dt_keep[choose][order], nu[choose][order], c=color, lw=1.0)
    ax.text(
            dt_keep[choose][0], nu[choose][0], "SN1998bw", fontsize=14,
            verticalalignment='center', horizontalalignment='right')


def sn2003L(ax, color):
    """ From Soderberg et al. 2005 """
    p = 3.2
    alpha_gamma = 0.075
    alpha_B = -1.0
    alpha_r = 0.96
    exp = (2*(p-2)*alpha_gamma + 2*(3+p/2)*alpha_B + 2*alpha_r)/(p+4)
    dt = np.linspace(44.5, 573.0)
    t0 = 30
    nu_a0 = 22.5
    nu_a = nu_a0 * (dt/t0)**(exp)
    ax.plot(dt, nu_a, c=color, lw=1.0)
    ax.text(
            dt[0], nu_a[0], "SN2003L", 
            horizontalalignment='right', fontsize=14,
            verticalalignment='center')
     

def sn2003bg(ax, color, label):
    """ From Soderberg et al. 2006 """
    p = 3.2
    alpha_gamma = -0.4
    alpha_B = -1.0
    alpha_r = 0.8
    exp = (2*(p-2)*alpha_gamma + 2*(3+p/2)*alpha_B + 2*alpha_r)/(p+4)
    dt = np.linspace(35.0, 351)
    t0 = 35
    nu_a0 = 22.5
    nu_a = nu_a0 * (dt/t0)**(exp)
    plt.scatter(dt[0], nu_a[0], c=color, marker='o')
    ax.plot(dt, nu_a, c=color, lw=1.0, label=label)
    ax.text(
            dt[0], nu_a[0], "SN2003bg", 
            horizontalalignment='left', fontsize=14,
            verticalalignment='bottom')


def swiftj1644(ax, color, label=None):
    """ 
    Plot for Swift J1644 
    Values directly from the paper
    """
    dt = [5, 10, 15, 22]
    nu_peak = [600, 250, 140, 80]
    ax.scatter(dt, nu_peak, marker='s', color=color)
    if label:
        ax.plot(dt, nu_peak, color=color, lw=1.0, label=label)
    else:
        ax.plot(dt, nu_peak, color=color, lw=1.0)
    ax.text(dt[0], nu_peak[0], "SwiftJ1644", fontsize=14)


def grb130427A(ax, color, label):
    """
    From Perley et al
    """
    dt = [10]
    nu_peak = [1E9]
    ax.scatter(dt, nu_peak, marker='s', color=color)
    if label:
        ax.plot(dt, nu_peak, color=color, lw=1.0, label=label)
    else:
        ax.plot(dt, nu_peak, color=color, lw=1.0)
    ax.text(dt[0], nu_peak[0], "GRB 130427A", fontsize=14)



def extra():
    """ other stuff """
    print("hi")
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


if __name__=="__main__":
    fig, ax = plt.subplots(1, 1, figsize=(8,6))

    swiftj1644(ax, '#a6cee3', 'Rel. TDE')

    at2018cow(ax)
    sn2009bb(ax, '#1f78b4', "Rel. SN")
    sn1998bw(ax, '#1f78b4', "Rel. SN")
    sn2003L(ax, '#33a02c')
    sn2003bg(ax, '#33a02c', 'Interacting SN')
    sn1993J(ax, '#33a02c')
    grb130427A(ax, 'black', 'GRB')

    ax.set_ylabel("$\\nu_p$ [GHz]", fontsize=16)
    ax.set_xlabel("$\\Delta t$ [days]", fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2, 100)
    ax.set_ylim(0.5, 1000)
    plt.tight_layout()
    plt.legend(fontsize=14)

    plt.show()
    #plt.savefig("nupeak_evolution.png")
