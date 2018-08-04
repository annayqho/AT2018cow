""" Plot nu_peak over time for various transients, including AT2018cow """


import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.table import Table
from plot_spec import fit_self_abs


def sn2009bb(ax):
    """ Plot for SN 2009bb """
    dt = [20, 52, 81, 145]
    nu_peak = [6, 3, 1, 0.8]
    ax.scatter(dt, nu_peak, c='grey', marker='o')
    ax.plot(dt, nu_peak, c='grey', lw=2.0)
    ax.text(
            dt[-1], nu_peak[-1], "SN2009bb", 
            horizontalalignment='left', fontsize=14,
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
    ax.text(
            dt[choose][-1], nu[choose][-1], "AT2018cow", fontsize=14,
            verticalalignment='center')


def get_peak(ax, freq_raw, flux_raw):
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
            first = freq[0:ind]
            second = freq[ind:]
            #print(freq[0:ind])
            b = fit_self_abs(freq[0:ind], flux[0:ind])
            yfit_up = 10**(2*np.log10(xgrid)+b)
            #ax.plot(xgrid, yfit_up, ls=':', c='k')

            # plot and fit a nu*something line
            m,b = np.polyfit(np.log10(freq[ind:]), np.log10(flux[ind:]), deg=1)
            yfit_down = 10**(m*np.log10(xgrid)+b)
            #ax.plot(xgrid, yfit_down, ls=':', c='k')

            nupeak = xgrid[np.argmin(np.abs(yfit_up-yfit_down))]
            return (nupeak, False)
 


def sn1998bw(ax):
    """ Plot for SN 1998bw 
    http://adsabs.harvard.edu/abs/1999A%26AS..138..467W
    """
    dat = Table.read(
            "../data/radio_compilations/1998bw.dat", 
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
            nupeak_val,islim_val = get_peak(ax, freq[choose], flux[choose])
            nu.append(nupeak_val)
            islim.append(islim_val)

    dt_keep = np.array(dt_keep)
    nu = np.array(nu)
    islim = np.array(islim)

    #choose = islim == True
    #ax.scatter(dt_keep[choose], nu[choose], c='grey', marker='^')
    #plt.show()

    choose = islim == False
    ax.scatter(dt_keep[choose], nu[choose], c='grey', marker='s')

    #ax.plot(dt_keep, nu, c='grey', lw=2.0)
    ax.text(
            dt_keep[choose][-1], nu[choose][-1], "SN1998bw", fontsize=14,
            verticalalignment='center')



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

    # 2003bg
    plt.scatter([23, 351], [25, 2.2], c='grey', marker='s')
    plt.plot([23, 351], [25, 2.2], c='grey', marker='s')
    plt.text(351, 2.2, "2003bg", fontsize=14)


    # 2003L
    # http://iopscience.iop.org/article/10.1086/427649/pdf
    plt.scatter([30], [22.5], c='grey', marker='s')
    plt.text(30, 22.5, "2003L")


if __name__=="__main__":
    fig, ax = plt.subplots(1, 1, figsize=(8,6))

    at2018cow(ax)
    sn2009bb(ax)
    sn1998bw(ax)

    # ax.set_ylabel("$\\nu_p$ [GHz]", fontsize=16)
    # ax.set_xlabel("$\\Delta t$ [days]", fontsize=16)
    # ax.tick_params(axis='both', labelsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_xlim(2, 2000)
    # ax.set_ylim(0.5, 400)
    # plt.tight_layout()

    plt.show()
    #plt.savefig("nupeak_evolution.png")