""" Plots of highly zoomed-in light curves during long tracks """

import matplotlib.pyplot as plt
from matplotlib import ticker
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
import numpy as np


def update(dt, nu, f, ef, dt_val, nu_val, f_val, ef_val):
    dt.append(dt_val)
    nu.append(nu_val)
    fmean, wsum = np.average(f_val, weights=1/ef_val**2, returned=True)
    efmean = np.sqrt(1/wsum)
    f.append(fmean)
    ef.append(efmean)
    return dt, nu, f, ef


def run_day(ax, dt, nu, f, ef, xticks, yticks, day, legpos='lower right'):
    choose = np.logical_or(nu == 215.5, nu == 243.3)
    col = '#440154' # dark purple

    if sum(choose) > 0:
        ax.errorbar(
                dt[choose], f[choose], ef[choose], 
                fmt='s', c=col, lw=1.5, 
                label="%s GHz" %nu[choose][0])
        ax.plot(
                dt[choose], f[choose], linestyle='-', c=col, lw=1.5) 

    choose = np.logical_or(nu == 231.5, nu == 259.3)
    col = '#5ec962' # light green
    if sum(choose) > 0:
        ax.errorbar(
                dt[choose], f[choose], ef[choose],
                fmt='*', c=col, lw=0.5, alpha=0.8, ms=13,
                label="%s GHz" %nu[choose][0])
        ax.plot(
                dt[choose], f[choose], linestyle='-', c=col, lw=1.5)

    choose = np.logical_or.reduce((nu == 330.8, nu == 344.8, nu == 341.5))
    col = '#3b528b' # dark blue
    if sum(choose) > 0:
        print(len(ef), len(choose))
        ax.errorbar(
                dt[choose], f[choose], ef[choose],
                fmt='o', c=col, lw=0.5, alpha=0.8, ms=7,
                label="%s GHz" %nu[choose][0])
        ax.plot(
                dt[choose], f[choose], linestyle='-', c=col, lw=1.5)

    choose = np.logical_or.reduce((nu == 357.5, nu==346.8, nu == 360.8))
    col = '#fde725' # yellow
    if sum(choose) > 0:
        ax.errorbar(
                dt[choose], f[choose], ef[choose],
                fmt='o', c=col, lw=0.5, alpha=0.8, ms=7,
                label="%s GHz" %nu[choose][0])
        ax.plot(
                dt[choose], f[choose], linestyle='-', c=col, lw=1.5)

    ax.text(
            0.05, 0.95, day, fontsize=10, 
            transform=ax.transAxes, verticalalignment='top')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.tick_params(axis='both', labelsize=14)
    ax.legend(fontsize=10, loc=legpos, ncol=2)
    ax.set_yticks(yticks)
    ax.set_xticks(xticks)
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.minorticks_off()


def jun21(ax):
    dt = []
    nu = []
    f = []
    ef = []

    # June 21: 215.5, 231.5
    # let's measure the times from UT 0
    dt, nu, f, ef = update(dt, nu, f, ef, 5, 215.5, 
            np.array([11.31, 7.34]), np.array([1.94, 2.39]))
    dt, nu, f, ef = update(dt, nu, f, ef, 5, 231.5, 
            np.array([13.07, 11.81]), np.array([2.33, 2.63]))
    dt, nu, f, ef = update(dt, nu, f, ef, (6+8.5)/2, 215.5, 
            np.array([10.48, 12.72]), np.array([1.13, 1.57]))
    dt, nu, f, ef = update(dt, nu, f, ef, (6+8.5)/2, 231.5, 
            np.array([12.41, 15.01]), np.array([1.33, 1.70]))
    dt, nu, f, ef = update(dt, nu, f, ef, (9+11.5)/2, 215.5, 
            np.array([15.44, 14.05]), np.array([1.11, 1.35]))
    dt, nu, f, ef = update(dt, nu, f, ef, (9+11.5)/2, 231.5, 
            np.array([14.96, 14.93]), np.array([1.31, 1.45]))
    dt, nu, f, ef = update(dt, nu, f, ef, (11.5+14)/2, 215.5, 
            np.array([15.28, 15.36]), np.array([1.27, 1.52]))
    dt, nu, f, ef = update(dt, nu, f, ef, (11.5+14)/2, 231.5, 
            np.array([13.43, 15.95]), np.array([1.53, 1.64]))

    dt = np.array(dt)
    nu = np.array(nu)
    f = np.array(f)
    ef = np.array(ef)

    yticks = [8, 10, 12, 14, 16]
    xticks = [4, 6, 8, 10, 12]
    run_day(ax, dt, nu, f, ef, xticks, yticks, 'Day 5', 'lower right')
    ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)


def jun22(ax):
    dt = []
    nu = []
    f = []
    ef = []

    # June 21: 215.5, 231.5
    # let's measure the times from UT 0
    dt, nu, f, ef = update(dt, nu, f, ef, 
            (6.5+7.4)/2, 215.5, 
            np.array([33.67, 35.73]), np.array([1.92, 2.29]))
    dt, nu, f, ef = update(dt, nu, f, ef, 
            (6.5+7.4)/2, 231.5, 
            np.array([36.83, 33.36]), np.array([2.18, 2.42]))
    dt, nu, f, ef = update(dt, nu, f, ef,
            (7.4+8.4)/2, 215.5, 
            np.array([38.12, 39.92]), np.array([1.84, 2.21]))
    dt, nu, f, ef = update(dt, nu, f, ef,
            (7.4+8.4)/2, 231.5, 
            np.array([35.79, 38.77]), np.array([2.10, 2.34]))

    dt = np.array(dt)
    nu = np.array(nu)
    f = np.array(f)
    ef = np.array(ef)

    yticks = [32, 34, 36, 38, 40]
    xticks = [6, 7, 8, 9]

    run_day(ax, dt, nu, f, ef, xticks, yticks, 'Day 6', legpos='lower right')


def jun23(ax):
    # June 23: RxA was 215.5, 231.5 and RxB was 330.8, 346.8
    dt_single = np.array(
            [(3.5+5.5)/2, (5.5+7.5)/2, (7.5+9.5)/2, (9.5+11.5)/2, (11.5+14)/2])
    dt = np.copy(dt_single)
    nu = np.array([215.5]*5)
    f = np.array([35.34, 35.96, 35.82, 38.16, 39.76])
    ef = np.array([1.45, 1.12, 1.18, 1.02, 1.19])

    dt = np.append(dt, dt_single)
    nu = np.append(nu, np.array([231.5]*5))
    f = np.append(f, np.array([45.12, 40.64, 33.80, 38.90, 36.03]))
    ef = np.append(ef, np.array([1.72, 1.30, 1.36, 1.19, 1.39]))

    dt = np.append(dt, dt_single)
    nu = np.append(nu, np.array([330.8]*5))
    f = np.append(f, np.array([36.91, 28.07, 44.04, 39.59, 26.73]))
    ef = np.append(ef, np.array([13.7, 4.72, 4.22, 3.96, 6.05]))

    dt = np.append(dt, dt_single)
    nu = np.append(nu, np.array([346.8]*5))
    f = np.append(f, np.array([10.39, 25.14, 32.71, 32.94, 32.44]))
    ef = np.append(ef, np.array([13.0, 4.45, 3.87, 3.42, 4.79]))

    yticks = [20, 26, 32, 38, 44, 50]
    xticks = [4, 6, 8, 10, 12, 14]
    run_day(ax, dt, nu, f, ef, xticks, yticks, 'Day 7', legpos='lower right')
    ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)


def jun24(ax):
    # June 24: RxB was (344.8, 360.8) and RxA was (215.5, 231.5)
    dt_single = np.array([(4+6)/2, (6+8)/2, (8+10)/2, (10+12)/2, (12+14)/2])
    dt = np.copy(dt_single)
    nu = np.array([215.5]*5)
    f = np.array([43.62, 40.10, 36.20, 41.92, 38.61])
    ef = np.array([1.10, 0.87, 1.00, 1.09, 1.18])

    dt = np.append(dt, dt_single) 
    nu = np.append(nu, np.array([231.5]*5))
    f = np.append(f, np.array([42.58, 39.69, 36.36, 42.39, 40.45]))
    ef = np.append(ef, np.array([1.29, 1.01, 1.18, 1.29, 1.39]))

    dt = np.append(dt, dt_single) 
    nu = np.append(nu, np.array([344.8]*5))
    f = np.append(f, np.array([27.15, 26.71, 21.39, 26.18, 26.54]))
    ef = np.append(ef, np.array([4.24, 2.70, 2.79, 2.53, 3.75]))

    dt = np.append(dt, dt_single) 
    nu = np.append(nu, np.array([360.8]*5))
    f = np.append(f, np.array([28.31, 21.07, 19.07, 20.91, 27.73]))
    ef = np.append(ef, np.array([5.34, 3.13, 3.15, 2.90, 4.63]))

    dt = np.array(dt)
    nu = np.array(nu)
    f = np.array(f)
    ef = np.array(ef)

    yticks = [10, 20, 30, 40, 50, 60]
    xticks = [4, 6, 8, 10, 12, 14]

    run_day(ax, dt, nu, f, ef, xticks, yticks, 'Day 8', legpos='lower left')
    ax.set_xlabel("Hours since 0h UT", fontsize=16)


def jun25(ax):
    # June 25: RxA were 341.5, 357.5, and RxB were 243.3, 259.3
    dt_single = np.array([(3.5+5.5)/2, (5.5+7.5)/2, (7.5+9.5)/2])
    dt = np.copy(dt_single)
    nu = np.array([341.5]*3)
    f = np.array([23.06, 26.12, 19.75])
    ef = np.array([3.96, 2.26, 2.57])

    dt = np.append(dt, dt_single)
    nu = np.append(nu, np.array([357.5]*3))
    f = np.append(f, np.array([27.12, 31.49, 22.54]))
    ef = np.append(ef, np.array([7.11, 3.83, 4.31]))

    dt = np.append(dt, dt_single)
    nu = np.append(nu, np.array([243.3]*3))
    f = np.append(f, np.array([36.74, 37.20, 38.20]))
    ef = np.append(ef, np.array([1.62, 1.24, 1.38]))

    dt = np.append(dt, dt_single)
    nu = np.append(nu, np.array([259.3]*3))
    f = np.append(f, np.array([38.10, 39.43, 33.62]))
    ef = np.append(ef, np.array([2.12, 1.59, 1.77]))

    xticks = [3, 5, 7, 9, 11]
    yticks = [10, 20, 30, 40, 50]
    run_day(ax, dt, nu, f, ef, xticks, yticks, 'Day 9', legpos='lower right')
    ax.set_xlabel("Hours since 0h UT", fontsize=16)
    ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)


if __name__=="__main__":
    fig,axarr = plt.subplots(3,2,figsize=(8,8))
    jun21(axarr[0,0])
    jun22(axarr[0,1])
    jun23(axarr[1,0])
    jun24(axarr[1,1])
    jun25(axarr[2,0])
    plt.subplots_adjust(hspace=0.2)
    axarr[2,1].axis('off')
    #plt.tight_layout()
    #plt.show()
    plt.savefig("lc_zoom.png")
