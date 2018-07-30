""" Peak lum phase space from Soderberg (2009) with AT2018cow
and Swift J1644 on top """

from matplotlib import rc
rc("font", family="serif")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

direc = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
inputf = direc + "/Soderberg2009Fig4.1"

fig = plt.figure(figsize=(7,6))

dat = Table.read(inputf, format='ascii')
x = dat['col1']
y = dat['col2']

choose = np.logical_and(y < 1e50, x < 0.8)
plt.scatter(
        x[choose], y[choose], marker='o', c='k', s=150,
        label="SNe Ibc")

choose = np.logical_and(x > 0.5, x < 1.2)
plt.scatter(
        x[choose], y[choose], marker='s', c='k', s=150,
        label="Rel. SNe")

choose = y > 7e49
plt.scatter(
        x[choose], y[choose], marker='o', edgecolors='k', s=150,
        facecolors='none', label="GRBs")


# AT2018cow
plt.scatter(
        0.173, 1.2E49, 
        marker='*', s=300, facecolors='white', edgecolors='black')
plt.text(
        0.173, 2.2E49, "AT2018cow", 
        fontsize=14, verticalalignment='bottom', horizontalalignment='center')


plt.legend()

plt.xlim(0.02, 12)
plt.ylim(1E45, 1E53)
plt.xscale('log')
plt.yscale('log')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel(
    "Blastwave Velocity $(\\Gamma \\beta)$",
    fontsize=14)
plt.ylabel(
    "Energy (erg)", fontsize=14)
plt.tight_layout()

#plt.show()
plt.savefig("soderberg_plot.png")
