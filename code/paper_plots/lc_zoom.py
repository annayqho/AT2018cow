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

    # PLOT
    dt = np.array(dt)
    nu = np.array(nu)
    f = np.array(f)
    ef = np.array(ef)

    choose = nu == 215.5
    ax.errorbar(
            dt[choose], f[choose], ef[choose], 
            fmt='s', c='k', lw=1.5, label="215.5 GHz")
    ax.plot(
            dt[choose], f[choose], linestyle='-', c='k', lw=1.5) 

    choose = nu == 231.5
    ax.errorbar(
            dt[choose], f[choose], ef[choose],
            fmt='*', c='#66c2a5', lw=0.5, alpha=0.8, ms=13,
            label="231.5 GHz")
    ax.plot(
            dt[choose], f[choose], linestyle='-', c='#66c2a5', lw=1.5)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)
    ax.set_xlabel("Hours since June 21, 0h UT", fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.legend(fontsize=12, loc='lower left')
    ax.set_yticks([8, 10, 12, 14, 16])
    ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks([4, 6, 8, 10, 12])
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.minorticks_off()


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

    # PLOT
    dt = np.array(dt)
    nu = np.array(nu)
    f = np.array(f)
    ef = np.array(ef)

    choose = nu == 215.5
    ax.errorbar(
            dt[choose], f[choose], ef[choose], 
            fmt='s', c='k', lw=1.5, label="215.5 GHz")
    ax.plot(
            dt[choose], f[choose], linestyle='-', c='k', lw=1.5) 

    choose = nu == 231.5
    ax.errorbar(
            dt[choose], f[choose], ef[choose],
            fmt='*', c='#66c2a5', lw=0.5, alpha=0.8, ms=13,
            label="231.5 GHz")
    ax.plot(
            dt[choose], f[choose], linestyle='-', c='#66c2a5', lw=1.5)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)
    ax.set_xlabel("Hours since June 22, 0h UT", fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.legend(fontsize=12, loc='lower left')
    ax.set_yticks([32, 34, 36, 38, 40])
    ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks([6, 7, 8, 9])
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.minorticks_off()


def jun23(ax):
    dt = []
    nu = []
    f = []
    ef = []

    # June 23: 215.5, 231.5
    # let's measure the times from UT 0
    dt, nu, f, ef = update(dt, nu, f, ef, 
            (3.5+5.5)/2, 215.5, 
            np.array([35.34, 36.91]), np.array([1.45, 13.7]))
    dt, nu, f, ef = update(dt, nu, f, ef, 
            (3.5+5.5)/2, 231.5, 
            np.array([45.12, 10.39]), np.array([1.72, 13.0]))
    dt, nu, f, ef = update(dt, nu, f, ef,
            (5.5+7.5)/2, 215.5, 
            np.array([35.96, 28.07]), np.array([1.12, 4.72]))
    dt, nu, f, ef = update(dt, nu, f, ef,
            (5.5+7.5)/2, 231.5, 
            np.array([40.63, 25.14]), np.array([1.30, 4.45]))
    dt, nu, f, ef = update(dt, nu, f, ef,
            (7.5+9.5)/2, 215.5, 
            np.array([35.82, 44.04]), np.array([1.18, 4.22]))
    dt, nu, f, ef = update(dt, nu, f, ef,
            (7.5+9.5)/2, 231.5, 
            np.array([33.80, 32.71]), np.array([1.36, 3.87]))
    dt, nu, f, ef = update(dt, nu, f, ef,
            (9.5+11.5)/2, 215.5, 
            np.array([38.16, 39.59]), np.array([1.02, 3.96]))
    dt, nu, f, ef = update(dt, nu, f, ef,
            (9.5+11.5)/2, 231.5, 
            np.array([38.90, 32.94]), np.array([1.19, 3.42]))

    # PLOT
    dt = np.array(dt)
    nu = np.array(nu)
    f = np.array(f)
    ef = np.array(ef)

    choose = nu == 215.5
    ax.errorbar(
            dt[choose], f[choose], ef[choose], 
            fmt='s', c='k', lw=1.5, label="215.5 GHz")
    ax.plot(
            dt[choose], f[choose], linestyle='-', c='k', lw=1.5) 

    choose = nu == 231.5
    ax.errorbar(
            dt[choose], f[choose], ef[choose],
            fmt='*', c='#66c2a5', lw=0.5, alpha=0.8, ms=13,
            label="231.5 GHz")
    ax.plot(
            dt[choose], f[choose], linestyle='-', c='#66c2a5', lw=1.5)

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel("$f_{\\nu}$ [mJy]", fontsize=16)
    ax.set_xlabel("Hours since June 23, 0h UT", fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.legend(fontsize=12, loc='lower left')
    ax.set_yticks([28, 32, 36, 40, 44, 48])
    ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks([4, 6, 8, 10, 12])
    ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax.minorticks_off()


if __name__=="__main__":
    fig,axarr = plt.subplots(2,2)
    jun21(axarr[0,0])
    jun22(axarr[0,1])
    jun23(axarr[1,0])
    plt.tight_layout()
    plt.show()
