""" Plot the X-ray light curve """

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
#from plot_spec import get_data
#from plot_lc import sma
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)

direc = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"

def get_xrt():
    """ Get the X-ray data
    X ray: 2018 Jun 19 at 10:34:30.742 UT was T0
    First SMA track was 2018 Jun 21, 03:58:25 - 14:18:09 UTC (10h 19m)
    I call this Day 5. That means that the XRT point is Day 3.
    So I need to add 3 days to all XRT points.
    """
    dat = Table.read(direc + "/xray_lc.txt", format='ascii')
    t = dat['col1'] / (3600*24) # in days
    # count to flux conversion (absorbed): 4.26 × 10-11 erg cm-2 ct-1
    flux = dat['col4'] * 4.26 * 10**(-11)
    eflux = dat['col5']
    t_xray = (t-t[0])+3
    return t_xray, flux, eflux


def get_nustar():
    """ Get the NuSTAR data provided by Varun
    Explosion date MJD 58285
    """
    dat = Table.read(direc + "/nustar.txt", format='ascii')
    dt = dat['mjd'] - 58285
    flux380 = 10**(dat['flux380'])
    flux380_low = 10**dat['flux380_low']
    flux380_high = 10**dat['flux380_high']
    flux380_err = (flux380_high-flux380_low)/2
    flux320 = 10**dat['flux320']
    flux320_low = 10**dat['flux320_low']
    flux320_high = 10**dat['flux320_high']
    flux320_err = (flux320_high-flux320_low)/2
    return dt, flux380, flux380_err, flux320, flux320_err



#def plot():
#     fig, axarr = plt.subplots(3, 1, sharex=True)
#     ax = axarr[0]
#     dat = Table.read(
#         "../data/radio_lc.dat", delimiter="&", format='ascii.no_header')
#     tel = np.array(dat['col2'])
#     choose = tel == 'SMA'
# 
#     days = np.array(dat['col1'][choose])
#     freq = np.array(dat['col3'][choose]).astype(float)
#     flux_raw = np.array(dat['col4'][choose])
#     flux = np.array(
#             [float(val.split("pm")[0][1:]) for val in flux_raw])
#     eflux = np.array(
#             [float(val.split("pm")[1][0:-1]) for val in flux_raw])
# 
#     choose = np.logical_and(freq > 230, freq < 245)
#     low_freq = min(freq[choose])
#     max_freq = max(freq[choose])
# 
#     f_radio = flux[choose] * 1e-3 * 1e-23 * freq[choose] * 1e9 / 1e-13
#     ef_radio = eflux[choose] * 1e-3 * 1e-23 * freq[choose] * 1e9 / 1e-13
#     t_radio = days[choose]
#     ax.errorbar(
#             t_radio, f_radio, yerr=ef_radio,
#             fmt='.', c='k', lw=0.5)
#     ax.plot(
#             t_radio, f_radio, linestyle='-', c='k',
#             label="%s-%s GHz" %(low_freq, max_freq), lw=1.0)
# 
#     ax.set_ylabel("$F_R = \\nu f_{\\nu}$", fontsize=16)
#     ax.yaxis.set_tick_params(labelsize=14)
# 
#     ax.text(0.05, 0.9, "SMA", transform=ax.transAxes,
#             fontsize=14, verticalalignment='top')
#     ax.text(0.9, 0.9, 
#     "$10^{-13} \mathrm{erg}\,\mathrm{s}^{-1}\,\mathrm{cm}^{-2}$", 
#     fontsize=14, transform=ax.transAxes, 
#     verticalalignment='top', horizontalalignment='right')
# 
# 
#     ax = axarr[1]
#     ax.errorbar(t_xray, f_xray, yerr=ects, fmt='.', c='k')
#     #ax.plot((t-t[0])+3, cts, c='grey', lw=2.0)
#     ax.set_ylabel(
#     r"$F_X = \nu f_\nu$", fontsize=16)
# 
#     ax.text(0.9, 0.9, 
#             "$10^{-11} \, \mathrm{erg \, s^{-1}} \, \mathrm{cm}^{-2}$",
#             fontsize=14, transform=ax.transAxes, 
#             verticalalignment='top', horizontalalignment='right')
#     ax.yaxis.set_tick_params(labelsize=14)
#     ax.xaxis.set_tick_params(labelsize=14)
# 
#     #for ii in [12.9, 21.4, 24.4]:
#     #    ax.axvline(x=ii, alpha=0.5, lw=2)
# 
#     ax.text(
#             0.1, 0.9, "Swift/XRT", 
#             transform=ax.transAxes, fontsize=14, verticalalignment='top')
# 
#     ax = axarr[2]
#     xgrid = np.linspace(5, 45)
#     ygrid = np.interp(xgrid, t_radio, f_radio*1e-13) / np.interp(xgrid, t_xray, f_xray*1e-11)
#     ax.plot(xgrid, ygrid, c='black')
#     ax.axvline(x=22, c='grey', ls='--')
#     ax.set_xlabel("Days Since 16 June UT", fontsize=14)
#     ax.set_ylabel(
#     r"$F_R/F_X$", fontsize=16)
#     ax.yaxis.set_tick_params(labelsize=14)
#     ax.xaxis.set_tick_params(labelsize=14)
# 
# 
#     #plt.show()
#     plt.savefig("xray_vs_radio.png")
