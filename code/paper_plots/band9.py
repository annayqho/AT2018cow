""" Light curve showing the weird Band 9 detection """

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from plot_spec import get_data

tel, freq, days, flux, flux_err = get_data()

fig = plt.figure(figsize=(4,4))

choose = np.logical_or(days == 22, days == 24)

plt.scatter(freq[choose], flux[choose], c='k')
plt.gca().tick_params(axis='both', labelsize=14)
plt.xlabel("Frequency [GHz]", fontsize=16)
plt.ylabel("Flux [mJy]", fontsize=16)
plt.xscale('log')
plt.yscale('log')

plt.show()
