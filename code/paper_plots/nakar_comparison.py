""" Compare the velocity-energy of the source to the shock breakout
parameters in Nakar & Sari """

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
import numpy as np

labels = ['RSG', 'BSG', 'WR'] # kind of star
v = [7000, 20000, 40000] # km/s
E = [9e47, 3e46, 9e44] # erg
colors = ['#d95f02', '#7570b3', '#1b9e77']

# AT2018cow
plt.scatter(
        0.14*3*10**5, 1.2e49, 
        marker='*', s=200, c='k', label='AT2018cow')

# SN 2003L
plt.scatter(
        (1/6)*3*10**5, 1.5e47, 
        marker='s', c='k', label='SN2003L')

# SN 2003bg
plt.scatter(
        0.13*3*10**5, 1.7e48, 
        marker='^', c='k', label='SN2003bg')

# SN 2007bg
plt.scatter(
        0.2*3*10**5, 1e48, 
        marker='v', c='k', label='SN2007bg')

for ii,lab in enumerate(labels):
    plt.scatter(v[ii], E[ii], c=colors[ii], label=lab)

plt.plot(v, E, c='k', lw=0.5, ls='--')


plt.legend(fontsize=14)
plt.xlabel("Velocity (km/s)", fontsize=16)
plt.ylabel("Energy in radio-emitting region (erg)", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim(1e44, 1e50)
plt.yscale('log')

plt.tight_layout()
plt.show()
#plt.savefig("vel_e_nakar.png")
