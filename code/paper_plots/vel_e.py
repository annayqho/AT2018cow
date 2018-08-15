""" Peak lum phase space from Soderberg (2009) with AT2018cow
and Swift J1644 on top """

from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

direc = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
inputf = direc + "/Soderberg2009Fig4.1"

fig = plt.figure(figsize=(6,5))

dat = Table.read(inputf, format='ascii')
x = dat['col1']
y = dat['col2']

choose = np.logical_and(y < 1e50, x < 0.8)
plt.scatter(
        x[choose], y[choose], marker='o', c='k', s=100,
        label="SNe Ibc")

choose = np.logical_and(x > 0.5, x < 1.2)
plt.scatter(
        x[choose], y[choose], marker='s', c='k', s=100,
        label="Rel. SNe")

choose = y > 7e49
plt.scatter(
        x[choose], y[choose], marker='o', edgecolors='k', s=100,
        facecolors='none', label="GRBs")


# AT2018cow
plt.scatter(
        0.173, 1.2E49, 
        marker='*', s=300, facecolors='white', edgecolors='black')
plt.text(
        0.173, 2.2E49, "AT2018cow", 
        fontsize=14, verticalalignment='bottom', horizontalalignment='center')


# Swift J1644
plt.scatter(
        0.5*1.2, 2.9E50 + 0.1*2.9E50, 
        marker='s', s=100, facecolors='white', edgecolors='black',
        label="Rel. TDE")

plt.legend(loc='upper left', fontsize=14)

plt.xlim(0.03, 8)
plt.ylim(1E45, 1E52)
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
