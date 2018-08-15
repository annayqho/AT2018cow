import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
import numpy as np


alma_bands = [3, 4, 5, 6, 7, 8, 9]
alma_freq = [98.9, 147, 206, 236, 348, 411, 689]

vla_bands = ['X', 'Ku', 'K', 'Ka', 'Q']
vla_freq = [10, 15, 22, 33, 45]

# plot these as light shaded regions
for ii,band in enumerate(vla_bands):
    plt.axvline(x=vla_freq[ii], color='#8da0cb', alpha=0.2, lw=7)
    plt.text(
             vla_freq[ii], 70, "VLA %s" %band, 
             rotation=90, horizontalalignment='center',
             verticalalignment='center') 

# plot the ATCA bands
# for ii,band in enumerate([9, 34]):
#     plt.axvline(x=band, color='#66c2a5', alpha=0.2, lw=7)
#     plt.text(
#             band, 80, "ATCA %sGHz" %band, 
#             rotation=90, horizontalalignment='center') 
# 

# plot the SMA bands
# for ii,band in enumerate([230, 340]):
#     plt.axvline(x=band, color='#fc8d62', alpha=0.2, lw=7)
#     plt.text(
#             band, 1.2, "SMA %sGHz" %band, 
#             rotation=90, horizontalalignment='center') 

freq = np.array([9, 34, 243, 259, 341, 357])
flux = np.array([0.46, 5.6, 35, 33, 17, 12])
flux_err = np.array([0.06, 0.16, 0.6, 0.75, 1.1, 2.0])

plt.errorbar(freq, flux, yerr=flux_err, fmt='.', c='k', marker='o')
plt.plot(freq, flux, linestyle='-', lw=0.5, c='k')

plt.text(
        freq[0]*1.1, flux[0], "(ATCA 9 GHz)", verticalalignment='top',
        horizontalalignment='left')

plt.text(
        freq[1]*1.1, flux[1], "(ATCA 34 GHz)", verticalalignment='top',
        horizontalalignment='left')

plt.text(
        freq[2]*1.1, flux[2], "(SMA 230 GHz)", verticalalignment='bottom',
        horizontalalignment='left')

plt.text(
        freq[5]*1.1, flux[5]*0.8, "(SMA 340 GHz)", verticalalignment='top',
        horizontalalignment='center')

# Plot nu**2
x = np.linspace(freq[0], freq[2])
y = (flux[0]-0.05) * (x/freq[0])**2 
plt.plot(x, y, linestyle='--', lw=1.0, c='k', alpha=0.2)
plt.text(12, 2.4, "$F_\\nu \propto \\nu^2$", fontsize=14)

# Plot nu**(1/3)
x = np.linspace(48, freq[2])
y = (24) * (x/70)**(1/3)
plt.plot(x, y, linestyle='--', lw=1.0, c='k', alpha=0.2)
plt.text(100, 33, "$F_\\nu \propto \\nu^{1/3}$", fontsize=14)

# Plot nu**(1/3)
#x = np.linspace(62, 100)
#y = 25 * (x/62)**(1/3)
#plt.plot(x, y, linestyle='--', lw=2.0, c='k')

#Plot nu**(-(p-1)/2)
# x = np.linspace(freq[2], 9e2)
# p=2
# y1 = flux[2] * (x/freq[2])**(-(p-1)/2)
# p=2.5
# y2 = flux[2] * (x/freq[2])**(-(p-1)/2)
# plt.fill_between(x, y1, y2, color='grey')
# 
# #Plot nu**(-p/2)
# x = np.linspace(freq[2], 9e2)
# p=2
# y1 = flux[2] * (x/freq[2])**(-p/2)
# p=2.5
# y2 = flux[2] * (x/freq[2])**(-p/2)
# plt.fill_between(x, y1, y2, color='grey')

# Plot the relativistic Maxwellian

#a = 9e3
#b = 0.2
#a = 5e2
#b = 0.05 
#x = np.linspace(1e5, 1e9)
#x = np.logspace(-1, 3)
#I = a * (2.5651*(1 + 1.92/((b*x)**(1/3)) + 0.9977/((b*x)**(2/3))) * np.exp(-1.8899*((b*x)**(1/3))))
#plt.plot(x,I,color='red')

plt.text(9, 170, "Spectrum at $\Delta t = 10\,$days", fontsize=14)

plt.xlabel("Frequency [GHz]", fontsize=16)
plt.ylabel("Flux [mJy]", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.xlim(7e0, 8e2)
plt.ylim(1e-1, 4e2)
plt.tight_layout()
#plt.show()

plt.savefig("vla_spec_day_10.png")
