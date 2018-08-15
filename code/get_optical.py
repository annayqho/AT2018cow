""" Get the optical data """

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

data_dir = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"

def get_opt_lum():
    dat = Table.read("%s/physevol.dat" %data_dir, format='ascii.no_header')
    mjd = dat['col1']
    lum = dat['col12'] * 4e33
    return mjd-mjd[0], lum


def get_opt_sed(day):
    filters = {'UVW2': 1928, 'UVM2': 2246, 'UVW1': 2600, 
               'u': 3543, 'g': 4770, 'r': 6231, 'i':7625, 'z': 9134,
               'H': 16350}
    dat = Table.read("%s/sed_uvoir.txt" %data_dir, format='ascii')
    choose = dat['Day']==day
    flux = []
    freq = []
    for col in dat.colnames[1:]:
        wl = filters[col]
        nu = (3e18 / wl)
        freq.append(nu)
        flux.append(dat[col][choose][0])
    return np.array(freq), np.array(flux)


if __name__=="__main__":
    print(get_opt_sed(5))
