""" Number of afterglows detectable n years later,
within say 200 Mpc, given the GRB (or LLGRB) rate """


import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
import numpy as np
from astropy.cosmology import Planck15


def get_years(redshift):
    """ return number of years source remains above threshold
    for a given distance """
    d = Planck15.luminosity_distance(z=redshift).cgs.value

    # initial flux is around 2e31 erg/s/Hz until 10-20 days post-burst
    # (Chandra & Frail 2012)
    f0 = 2e31 / (4 * np.pi * d**2)
    t_yr = np.logspace(-1, 4, 1000)
    f = f0 * ((t_yr*365) / (15))**(-1)

    # convert to mJy
    f_mjy = (f / 10**(-23)) * 10**3

    # let's say that "detectable" means above 100 uJy
    choose = f_mjy > 0.1

    return t_yr[choose][-1]


def make_grbs(t, redshift):
    """ Given the rate of GRBs out to a given redshift,
    estimate the number per time """


if __name__=="__main__":
    # rate of bursts 
    R_raw = 3e-8 # /yr/Mpc3

    z = np.linspace(0.01, 0.05)
    d_mpc = Planck15.luminosity_distance(z=z).value
    vol = np.array([(4/3) * np.pi * val**3 for val in d_mpc])
    R = R_raw * vol
    dur = np.array([get_years(val) for val in z])
    plt.plot(R, dur, c='k')
    #plt.xscale('log')
    #plt.yscale('log')
    plt.ylabel("Duration (years)", fontsize=14)
    plt.xlabel("Rate (N/yr)", fontsize=14)
    plt.show()


    ## there's around 1 burst per year in the local universe??
