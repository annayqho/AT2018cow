""" Plot the spectrum from radio to X-ray frequencies """

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
import numpy as np

# There is Swift/XRT data on June 19, June 21, June 29, and June 30
dates = ['June 21', 'June 29', 'June 30']
days = [5, 13, 14]
Fx_all = [1.2e-11, 1.8e-11, 6.05e-12]
fr_all = [13.26, 32, 43.5]
nu_x = 0.25e18
nu_r = 230e9
markers = ['o', 's', '^']
styles = ['--', '-.', ':']
try_index = 1.5

for ii,date in enumerate(dates):
    # Let's start with June 21:
    Fx = Fx_all[ii] # erg/cm2/s
    # 1 keV ~ 0.25e18 Hz
    fx = Fx/nu_x # erg/cm2/s/Hz
    fx_mJy = fx * 10**(23) * 10**3 # mJy 

    # On this day, our SMA flux at 230 GHz was 13.26 mJy
    fr_mJy = fr_all[ii]

    # Solve for the spectral index
    alpha = np.round(np.log10(fx_mJy/fr_mJy)/np.log10(nu_r/nu_x), 2)

    plt.scatter(nu_x, fx_mJy, c='k', marker=markers[ii])
    plt.scatter(nu_r, fr_mJy, c='k', marker=markers[ii])
    plt.plot(
            [nu_x, nu_r], [fx_mJy, fr_mJy], 
            ls=styles[ii], c='k', 
            label="Day %s: $F_\\nu \propto \\nu^{-%s}$" %(days[ii], alpha))
    if ii >= 1:
        plt.plot(
                [nu_r, nu_x],
                [fr_mJy, fr_mJy * (nu_r/nu_x)**try_index],
                c='k', ls='-', lw=0.5)
        plt.text(0.5e12, 0.1, "$F_\\nu \propto \\nu^{-1.5}$", fontsize=12)

plt.xscale('log')
plt.yscale('log')
plt.xlabel("Frequency [Hz]", fontsize=16)
plt.ylabel("Flux [mJy]", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()

plt.xlim(5e10, 1e18)
plt.ylim(0.0005, 100)
plt.tight_layout()

#plt.savefig("xray_spec_index.png")
plt.show()
