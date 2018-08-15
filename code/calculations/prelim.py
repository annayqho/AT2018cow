""" preliminary calculations """

from astropy.cosmology import Planck15
import numpy as np

# distance to the source in cm
d = Planck15.luminosity_distance(z=0.014).cgs.value
d27 = d/10**27

# Peak 1: 230 GHz, peak at 8 days at 40 mJy
F1 = 40
nu1 = 230

# Peak 2: 340 GHz, peak at 7 days at 36 mJy
F2 = 36
nu2 = 340

# Calculation
Fratio = F1/F2
nuratio = nu2/nu1
p = 1 + 2 * np.log10(Fratio) / np.log10(nuratio)
print(p)
