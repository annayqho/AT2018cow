""" 
Plot the spectrum from radio to X-ray frequencies 

Choose Days 10 (June 26), 13/14 (June 29/30), and 22 (July 8),
because these have the best sampling.
"""

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import numpy as np
from astropy.cosmology import Planck15
from xray_lc import get_xray
from plot_spec import get_data
from get_optical import get_opt_data


def get_sed(day):
    """ Return freq and nufnu for an individual day """
    freq = []
    nufnu = []

    # X rays
    freq.append(1E17)
    nufnu.append(np.interp(day, t_x, lum_x))
    
    # radio
    choose = t_rad == day
    for ii,val in enumerate(nu_rad[choose]):
        if ferr_rad[choose][ii] > 0:
            freq.append(val)
            nufnu.append(lum_rad[choose][ii]) 

    # optical
    # goes from 1,000 -- 20,000 A
    # which is 10^14 -- 10^15 Hz
    # I need to ask Dan for the SED
    freq.append(1E14)
    nufnu.append(np.interp(day, t_x, lum_x))
    return freq, nufnu



# set up the plot
fig = plt.figure(figsize=(6,4))

# load data
t_x, f_x, ef_x = get_xray()
d = Planck15.luminosity_distance(z=0.014).cgs.value
lum_x = f_x * 4 * np.pi * d**2
tel, nu_rad, t_rad, f_rad, ferr_rad = get_data()
nu_rad = nu_rad * 1E9
lum_rad = nu_rad * (f_rad*1e-3*1e-23) * 4 * np.pi * d**2
# need to get the optical LC
t_o, lum_o = get_opt_data()

# plot
x, y = get_sed(10)
plt.scatter(x, y, c='k', marker='o')
#x, y = get_sed(22) # and include Day 24 Band 9
#plt.scatter(x, y, c='grey', marker='s')
#x, y = get_sed(14)
#plt.scatter(x, y, c='grey', marker='o')


plt.xscale('log')
plt.yscale('log')
plt.xlabel("Frequency [Hz]", fontsize=16)
plt.ylabel(r"$\nu L_\nu$ (erg/s)", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()

plt.tight_layout()

#plt.savefig("xray_spec_index.png")
plt.show()
