""" Get the optical data """

import matplotlib.pyplot as plt
from astropy.table import Table


def get_opt_data():
    dat = Table.read("../data/physevol.dat", format='ascii.no_header')
    mjd = dat['col1']
    lum = dat['col12'] * 4e33
    return mjd-mjd[0], lum
