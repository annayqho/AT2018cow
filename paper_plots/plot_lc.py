""" Make a grid of light curves from the SMA, ATCA, and Swift/XRT """

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
import matplotlib.gridspec as gridspec
import numpy as np
from astropy.table import Table


def spindex(ax, d, nulow, nuhigh, flow, fhigh):
    alpha = np.log10( fhigh / flow ) / np.log10( nulow / nuhigh )
    ax.scatter(d, alpha, c='k')
    ax.plot(d, alpha, c='k', lw=1.0)
    ax.set_ylabel("$\\alpha$", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)


def sma(ax):
    dat = Table.read(
        "../data/radio_lc.dat", delimiter="&", format='ascii.no_header')
    tel = np.array(dat['col2'])
    choose = tel == 'SMA'

    days = np.array(dat['col1'][choose])
    freq = np.array(dat['col3'][choose]).astype(float)
    flux_raw = np.array(dat['col4'][choose])
    flux = np.array(
            [float(val.split("pm")[0][1:]) for val in flux_raw])
    flux_err = np.array(
            [float(val.split("pm")[1][0:-1]) for val in flux_raw])

    choose = np.logical_and(freq > 230, freq < 235)
    low_freq = min(freq[choose])
    max_freq = max(freq[choose])

    ax.errorbar(
            days[choose], flux[choose], yerr=flux_err[choose], 
            mfc='white', mec='black', 
            fmt='.', marker='o', label="%s-%s GHz" %(low_freq, max_freq), c='k')

    choose = np.logical_and(freq > 240, freq < 245)
    low_freq = min(freq[choose])
    max_freq = max(freq[choose])

    ax.errorbar(
            days[choose], flux[choose], yerr=flux_err[choose], 
            mfc='white', mec='black', 
            fmt='.', marker='s', label="%s-%s GHz" %(low_freq, max_freq), c='k')

    choose = np.logical_and(freq > 230, freq < 245)
    ax.plot(
            days[choose], flux[choose], linestyle='--', c='k')
    dlow_return = days[choose]
    flow_return_temp = flux[choose]
    nulow_return_temp = freq[choose]

    choose = np.logical_and(freq > 330, freq < 334)
    low_freq = min(freq[choose])
    max_freq = max(freq[choose])

    ax.errorbar(
            days[choose], flux[choose], yerr=flux_err[choose], 
            mfc='white', mec='black', ms=8,
            fmt='.', marker='v', label="%s-%s GHz" %(low_freq, max_freq), c='k')

    choose = np.logical_and(freq > 357, freq < 361)
    low_freq = min(freq[choose])
    max_freq = max(freq[choose])

    ax.errorbar(
            days[choose], flux[choose], yerr=flux_err[choose], 
            mfc='white', mec='black', ms=8,
            fmt='.', marker='^', label="%s-%s GHz" %(low_freq, max_freq), c='k')

    choose = np.logical_or(
            np.logical_and(freq > 330, freq < 334),
            np.logical_and(freq > 357, freq < 361))

    ax.plot(
            days[choose], flux[choose], linestyle='-', c='k')
    days_return = days[choose]
    fhigh_return = flux[choose]
    nuhigh_return = freq[choose]
    ind = np.array([np.where(dlow_return==val)[0][0] for val in days_return])
    flow_return = flow_return_temp[ind]
    nulow_return = nulow_return_temp[ind]

    ax.set_ylabel("$F_{\\nu}$ [mJy]", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)

    ax.legend()

    return days_return, nulow_return, nuhigh_return, flow_return, fhigh_return


def atca(ax):
    dat = Table.read(
        "../data/radio_lc.dat", delimiter="&", format='ascii.no_header')
    tel = np.array(dat['col2'])
    choose = tel == 'ATCA'

    days = np.array(dat['col1'][choose])
    freq = np.array(dat['col3'][choose]).astype(float)
    flux_raw = np.array(dat['col4'][choose])
    flux = np.array(
            [float(val.split("pm")[0][1:]) for val in flux_raw])
    flux_err = np.array(
            [float(val.split("pm")[1][0:-1]) for val in flux_raw])

    choose = freq==5.5

    ax.errorbar(
            days[choose], flux[choose], yerr=flux_err[choose], 
            mfc='white', mec='black', 
            fmt='.', marker='o', label="5.5 GHz", c='k')

    ax.plot(
            days[choose], flux[choose], linestyle='-', c='k')

    choose = freq == 9

    ax.errorbar(
            days[choose], flux[choose], yerr=flux_err[choose], 
            mfc='white', mec='black', 
            fmt='.', marker='s', label="9 GHz", c='k')

    ax.plot(
            days[choose], flux[choose], linestyle='--', c='k')


    choose = freq == 34

    ax.errorbar(
            days[choose], flux[choose], yerr=flux_err[choose], 
            mfc='white', mec='black', ms=8,
            fmt='.', marker='v', label="34 GHz", c='k')

    ax.plot(
            days[choose], flux[choose], linestyle=':', c='k')

    ax.set_xlabel("Days Since Discovery", fontsize=16)
    ax.set_ylabel("$F_{\\nu}$ [mJy]", fontsize=16)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)

    ax.set_yscale('log')

    ax.legend()


fig = plt.figure(figsize=(6,8))#, tight_layout=True)
gs = gridspec.GridSpec(2, 1, height_ratios=[3,1], hspace=0.1)
gs.update(left=0.15, right=0.9)

atca_ax = plt.subplot(gs[1])#fig.add_subplot(gs[5:6, :], sharex=ax1)
atca(atca_ax)

gs0 = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec = gs[0], height_ratios=[1,4], hspace=0)
#gs0.update(hspace=0.05)
ax1 = plt.subplot(gs0[1], sharex=atca_ax)#fig.add_subplot(gs[0:5, :])
days, nulow, nuhigh, flow, fhigh = sma(ax1)
plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = plt.subplot(gs0[0])#fig.add_subplot(gs[0:5, :])
spindex(ax2, days, nulow, nuhigh, flow, fhigh)
plt.setp(ax2.get_xticklabels(), visible=False)

plt.savefig("lc.png")
#plt.show()
