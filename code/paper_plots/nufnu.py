""" Estimate nu * f_nu for the radio data """

import numpy as np
import matplotlib.pyplot as plt

def get_flux(freq):
    peak_freq = 92 *1e9
    peak_flux = 92 *1e-3 * 1e-23
    flux = np.zeros(len(freq))
    choose = freq <= peak_freq
    flux[choose] = peak_flux * (freq[choose]/peak_freq)**2
    choose = freq > peak_freq
    flux[choose] = peak_flux * (freq[choose]/peak_freq)**(-1.2)
    return flux


# on Day 22

nu = np.linspace(10e9, 1000e9)
flux = get_flux(nu)

plt.plot(nu, nu*flux)
plt.xlim(10e9, 1000e9)
plt.show()


