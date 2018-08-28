""" Analysis of the spectral index of AT2018cow """

import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
import emcee
import corner
from get_radio import get_data_all, get_spectrum


tel, freq, day, flux, eflux_form, eflux_sys = get_data_all()

# Day 10
print("Analysis of Day 10")
choose_tel = tel == 'SMA'
choose_freq = freq > 100
choose_day = day == 10
choose = choose_tel & choose_freq & choose_day

nu = freq[choose]
print("Frequency")
print(nu)
f = flux[choose]
print("Flux")
print(f)
print("Unc. on Flux")
ef = np.sqrt(eflux_form[choose]**2 + eflux_sys[choose]**2)
print(ef)

npts = len(nu)


#plt.show()
plt.savefig("mcmc_fit_day_14.png")
