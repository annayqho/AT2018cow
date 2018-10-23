""" Get the optical data """

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.modeling.blackbody import blackbody_nu
from astropy.cosmology import Planck15

data_dir = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
d = Planck15.luminosity_distance(z=0.014).cgs.value


def get_opt_lum():
    dat = Table.read("%s/physevol.dat" %data_dir, format='ascii.no_header')
    mjd = dat['col1']
    lum = dat['col12'] * 4e33
    return mjd-mjd[0], lum


def get_opt_sed(day):
    filters = {'UVW2': 1928, 'UVM2': 2246, 'UVW1': 2600, 
               'u': 3543, 'g': 4770, 'r': 6231, 'i':7625, 'z': 9134,
               'H': 16350, 'K': 21900}
    dat = Table.read("%s/sed_uvoir.txt" %data_dir, format='ascii')
    choose = dat['Day']==day
    flux = []
    freq = []
    for col in dat.colnames[1:]:
        wl = filters[col]
        nu = (3e18 / wl)
        freq.append(nu)

        # ab mag
        ab = dat[col][choose][0]
        flux.append(10**(ab/(-2.5)) * 3631 * 1e3)
        # return in mJy
    return np.array(freq), np.array(flux)


def get_bb(day):
    # Get the blackbody curve on a certain day
    # what is the temp on that day?
    dat = Table.read("%s/physevol.dat" %data_dir, format='ascii.no_header')
    dt = dat['col2']
    Ts = dat['col6'] * 10**3
    T = Ts[np.argmin(np.abs(dt-day))]

    # what is the radius on that day?
    Rs = dat['col3'] * 1.5E13 # originally in AU
    R = Rs[np.argmin(np.abs(dt-day))]

    low_wl = 1928 # Angstroms, UVW2
    hi_wl = 21900 # Angstroms, K-band
    c = 3E18
    h = 4.136E-15
    k = 1.38E-16
    x = np.logspace(np.log10(c/hi_wl), np.log10(c/low_wl))
    y = blackbody_nu(x, T)
    print("temp", T)
    return x, y, R


def get_nonthermal(day):
    # Get the nonthermal component on a certain day
    # F_nu \propto (constant)*nu^{-0.75}
    # L2 is the nonthermal luminosity, which is column 15
    dat = Table.read("%s/physevol.dat" %data_dir, format='ascii.no_header')
    dt = dat['col2']
    Ls = dat['col15'] * 3.83E33 # originally in solar luminosities
    L = Ls[np.argmin(np.abs(dt-day))]
    F = L/(4*np.pi*d**2)

    # Now normalize the spectrum to that value
    low_wl = 1928 # Angstroms, UVW2
    hi_wl = 21900 # Angstroms, K-band
    c = 3E18
    low_freq = c/hi_wl
    hi_freq = c/low_wl
    # normalization constant
    K = F * 0.25 / (hi_freq**(0.25)-low_freq**(0.25))

    x = np.logspace(np.log10(c/hi_wl), np.log10(c/low_wl))
    y = K * x**(-0.75)
    return x, y # returning F_nu


if __name__=="__main__":
    print(get_opt_sed(5))
