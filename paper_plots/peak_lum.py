""" Peak lum phase space from Soderberg (2009) with AT2018cow
and Swift J1644 on top """

from matplotlib import rc
rc("font", family="serif")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

direc = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
inputf = direc + "/Soderberg2009Fig3.1"

fig = plt.figure(figsize=(7,6))

dat = Table.read(inputf, format='ascii')
x = dat['col1']
y = dat['col2']

choose = np.logical_or(y < 9E27, x > 50)
plt.scatter(
        x[choose], y[choose], marker='o', c='k', s=150,
        label="SNe Ibc")

choose = np.logical_and(y > 8E27, x < 80)
plt.scatter(
        x[choose], y[choose], marker='s', c='k', s=150,
        label="GRBs/Rel. SNe")

# Lines
# v = 0.1c
plt.plot([4.3,991], [1E25,9E29], ls='--', c='k')
plt.text(
        8, 2E25, "$R/\Delta t = 0.1c$", fontsize=14, rotation=45,
        horizontalalignment='left', verticalalignment='bottom')

# v = 0.01c
plt.plot([43,1000], [1E25, 7.4E27], ls='--', c='k')
plt.text(
        80, 2E25, "$R/\Delta t = 0.01c$", fontsize=14, rotation=45,
        horizontalalignment='left', verticalalignment='bottom')

# v = c
plt.plot([1, 207], [7E25, 4.4E30], ls='--', c='k')
plt.text(
        2, 1.8E26, "$R/\Delta t = c$", fontsize=14, rotation=45,
        horizontalalignment='left', verticalalignment='bottom')

# AT2018cow
plt.scatter(
        22*100/5, 2E29, 
        marker='*', s=300, facecolors='white', edgecolors='black')
plt.text(
        22*100/5, 3E29, "AT2018cow", 
        fontsize=14, verticalalignment='bottom', horizontalalignment='center')

plt.legend()

plt.xlim(1, 1000)
plt.ylim(1E25, 5E30)
plt.xscale('log')
plt.yscale('log')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel(
    "$(\Delta t/1\,\mathrm{day})(\\nu_p/5\,\mathrm{GHz})$",
    fontsize=14)
plt.ylabel("Peak Radio Luminosity ($\mathrm{erg\,s^{-1}\,Hz^{-1}}$)",
    fontsize=14)
plt.tight_layout()

#plt.show()
plt.savefig("chevalier_plot.png")
