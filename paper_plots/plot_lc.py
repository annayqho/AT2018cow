import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
import matplotlib.gridspec as gridspec
import numpy as np
from astropy.table import Table


def spindex(ax, d, nulow, nuhigh, flow, eflow, fhigh, efhigh):
    alpha = np.log10( fhigh / flow ) / np.log10( nulow / nuhigh )

    # estimate the uncertainty
    a = (1 / np.log10(nulow / nuhigh)) 
    b = flow/fhigh
    c = np.log10(np.e)
    dadf_hi = a*b*c* (1/flow) * efhigh
    dadf_low = a*b*c* (-fhigh/flow**2) * eflow
    ealpha = np.sqrt(dadf_hi**2 + dadf_low**2)

    #ax.scatter(d, alpha, c='k')
    ax.errorbar(d, alpha, yerr=ealpha, c='k', fmt='.')
    ax.plot(d, alpha, c='k', lw=1.0)
    ax.set_ylabel("$\\alpha$", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)

    return alpha, ealpha


def sma(ax, legend):
    dat = Table.read(
        "../data/radio_lc.dat", delimiter="&", format='ascii.no_header')
    tel = np.array(dat['col2'])
    choose = tel == 'SMA'

    days = np.array(dat['col1'][choose])
    freq = np.array(dat['col3'][choose]).astype(float)
    flux_raw = np.array(dat['col4'][choose])
    flux = np.array(
            [float(val.split("pm")[0][1:]) for val in flux_raw])
    eflux = np.array(
            [float(val.split("pm")[1][0:-1]) for val in flux_raw])

    # add 10% systematic uncertainty to everything
    flux_err = np.sqrt(eflux**2 + (0.10 * flux)**2)

    choose = np.logical_and(freq > 230, freq < 245)
    low_freq = min(freq[choose])
    max_freq = max(freq[choose])

    ax.errorbar(
            days[choose], flux[choose], flux_err[choose], 
            fmt='.', c='k', lw=0.5)
    ax.plot(
            days[choose], flux[choose], linestyle='-', c='k',
            label="%s-%s GHz" %(low_freq, max_freq))
    dlow_return = days[choose]
    flow_return_temp = flux[choose]
    eflow_return_temp = flux_err[choose]
    nulow_return_temp = freq[choose]

    choose = np.logical_and(freq > 330, freq < 361)
    low_freq = min(freq[choose])
    max_freq = max(freq[choose])

    ax.errorbar(
            days[choose], flux[choose], flux_err[choose],
            fmt='.', 
            label="%s-%s GHz" %(low_freq, max_freq), c='k')
    #ax.plot(
    #        days[choose], flux[choose], linestyle='-', c='k')
    days_return = days[choose]
    fhigh_return = flux[choose]
    efhigh_return = flux_err[choose]
    nuhigh_return = freq[choose]
    ind = np.array([np.where(dlow_return==val)[0][0] for val in days_return])
    flow_return = flow_return_temp[ind]
    eflow_return = eflow_return_temp[ind]
    nulow_return = nulow_return_temp[ind]

    ax.set_ylabel("$F_{\\nu}$ [mJy]", fontsize=16)
    ax.yaxis.set_tick_params(labelsize=14)

    ax.text(0.05, 0.9, "SMA", transform=ax.transAxes, fontsize=14)

    if legend:
        ax.legend()

    return days_return, nulow_return, nuhigh_return, flow_return, eflow_return, fhigh_return, efhigh_return


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
    eflux = np.array(
            [float(val.split("pm")[1][0:-1]) for val in flux_raw])

    # add 10% systematic uncertainty to everything
    flux_err = np.sqrt(eflux**2 + (0.10 * flux)**2)
    flux_err[eflux < 0] = -99 # non-detections

    choose_freq = freq==5.5

    # for detections
    det = flux_err > 0
    choose = np.logical_and(choose_freq, det)

    ax.errorbar(
            days[choose], flux[choose], flux_err[choose],
            mfc='white', mec='black', 
            fmt='.', marker='o', label="5.5 GHz", c='k')

    # for non-detections
    choose = np.logical_and(choose_freq, ~det)

    ax.scatter(
            days[choose], flux[choose], c='k', marker='v')

    # plot a line between all of them to make it clear
    # which frequency it corresponds to
    ax.plot(
            days[choose_freq], flux[choose_freq], linestyle='-', c='k')

    choose = freq == 9

    ax.errorbar(
            days[choose], flux[choose], flux_err[choose],
            mfc='white', mec='black', 
            fmt='.', marker='s', label="9 GHz", c='k')

    ax.plot(
            days[choose], flux[choose], linestyle='--', c='k')

    choose = freq == 34

    ax.errorbar(
            days[choose], flux[choose], flux_err[choose],
            mfc='white', mec='black', ms=8,
            fmt='.', marker='v', label="34 GHz", c='k')

    ax.plot(
            days[choose], flux[choose], linestyle=':', c='k')

    ax.set_xlabel("Days Since Discovery", fontsize=16)
    ax.set_ylabel("$F_{\\nu}$ [mJy]", fontsize=16)
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)

    ax.set_yscale('log')
    ax.text(
            0.05, 0.9, "ATCA", transform=ax.transAxes, fontsize=14,
            verticalalignment='top')

    ax.legend()


if __name__=="__main__":
    fig = plt.figure(figsize=(8,8))#, tight_layout=True)
    gs = gridspec.GridSpec(2, 1, height_ratios=[3,1], hspace=0.1)
    gs.update(left=0.15, right=0.9)

    atca_ax = plt.subplot(gs[1])#fig.add_subplot(gs[5:6, :], sharex=ax1)
    atca(atca_ax)

    gs0 = gridspec.GridSpecFromSubplotSpec(
            3, 1, subplot_spec = gs[0], height_ratios=[1,3,1], hspace=0)
    #gs0.update(hspace=0.05)
    ax1 = plt.subplot(gs0[1], sharex=atca_ax)#fig.add_subplot(gs[0:5, :])
    days, nulow, nuhigh, flow, eflow, fhigh, efhigh = sma(ax1, legend=True)
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax2 = plt.subplot(gs0[0], sharex=atca_ax)#fig.add_subplot(gs[0:5, :])
    alpha, ealpha = spindex(ax2, days, nulow, nuhigh, flow, eflow, fhigh, efhigh)
    plt.setp(ax2.get_xticklabels(), visible=False)

    #plt.savefig("lc.png")
    plt.show()
