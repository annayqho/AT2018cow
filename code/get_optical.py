""" Get the optical data """

import matplotlib.pyplot as plt
from astropy.table import Table


def get_opt_data():
    data_dir = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
    dat = Table.read("%s/physevol.dat" %data_dir, format='ascii.no_header')
    mjd = dat['col1']
    lum = dat['col12'] * 4e33
    return mjd-mjd[0], lum
