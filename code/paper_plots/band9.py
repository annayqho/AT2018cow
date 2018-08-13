""" Light curve showing the weird Band 9 detection """

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/code")
from get_radio import get_data_all, get_spectrum


fig = plt.figure(figsize=(6,4))

# Plot the Band 9 point
tel, freq, days, flux, eflux_form, eflux_sys = get_data_all()
choose = freq > 600
plt.errorbar(
        freq[choose], flux[choose], yerr=eflux_sys[choose], 
        fmt='o', c='k')

nu,flux = get_spectrum(24)
plt.scatter(nu, flux, c='k', marker='.')
plt.gca().tick_params(axis='both', labelsize=14)
plt.xlabel("Frequency [GHz]", fontsize=16)
plt.ylabel("Flux [mJy]", fontsize=16)
plt.xscale('log')
plt.yscale('log')

plt.show()
