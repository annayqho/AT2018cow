""" Format the 1993J measurements """

import numpy as np
from astropy.table import Table


def run(dt, freq, flux, dat, col, freqval):
    dt_temp = dat['col2']
    flux_temp = dat[col] 
    npts = len(flux_temp)
    dt.extend(list(dt_temp))
    for i in np.arange(npts):
        freq.append(float(freqval))
    flux.extend(list(flux_temp))
    return dt, freq, flux


def format_1993J():
    dt = [] # days
    freq = [] # GHz
    flux = [] # mJy

    # all these values are in mJy
    data_dir = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
    dat = Table.read(
        "%s/radio_compilations/1993J.dat" %data_dir, 
        delimiter='&', guess=False, format='ascii.no_header')
    dt_temp = dat['col2']

    dt, freq, flux = run(dt, freq, flux, dat, 'col4', 1.5)
    run(dt, freq, flux, dat, 'col5', 5.0)
    run(dt, freq, flux, dat, 'col6', 8.3)
    run(dt, freq, flux, dat, 'col7', 15.0)
    run(dt, freq, flux, dat, 'col8', 25)

    # extra values, including those at high frequency
    dat = Table.read(
        "../data/radio_compilations/1993J_extra.txt", 
        delimiter='&', guess=False, format='ascii.no_header')

    dt.extend(list(dat['col2']))
    freq.extend(list(dat['col5']))
    flux.extend(list(dat['col4']))

    dt = np.array(dt)
    freq = np.array(freq, dtype=float)
    flux = np.array(flux)

    # now go through and format the flux array

    flux[flux=='nodata'] = -99
    choose = np.array(['<' in val for val in flux])
    flux[choose] = -99

    for ii,fluxval in enumerate(flux):
        if '\\pm' in fluxval:
            flux[ii] = fluxval.split('\\pm')[0][1:].strip()
        else:
            if 'pm' in fluxval:
                flux[ii] = fluxval.split('pm')[0][1:].strip()

    flux = flux.astype(float)

    return dt, freq, flux
