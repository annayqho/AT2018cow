import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
import numpy as np

freq = np.array([9, 34, 243, 259, 341, 357])
flux = np.array([0.46, 5.6, 35, 33, 17, 12])
flux_err = np.array([0.06, 0.16, 0.6, 0.75, 1.1, 2.0])

plt.errorbar(freq, flux, yerr=flux_err, fmt='.', c='k', marker='o')
plt.plot(freq, flux, linestyle='-', lw=0.5, c='k')

# Plot nu**2
x = np.linspace(freq[0], freq[2])
y = (flux[0]-0.05) * (x/freq[0])**2 
plt.plot(x, y, linestyle='--', lw=1.0, c='k', alpha=0.2)
plt.text(12, 2.4, "$F_\\nu \propto \\nu^2$", fontsize=14)

# Plot nu**(1/3)
x = np.linspace(48, freq[2])
y = (24) * (x/70)**(1/3)
plt.plot(x, y, linestyle='--', lw=1.0, c='k', alpha=0.2)
plt.text(113, 33, "$F_\\nu \propto \\nu^{1/3}$", fontsize=14)

# Plot nu**(1/3)
#x = np.linspace(62, 100)
#y = 25 * (x/62)**(1/3)
#plt.plot(x, y, linestyle='--', lw=2.0, c='k')

#Plot nu**(-(p-1)/2)
x = np.linspace(freq[2], freq[-1])
p=2
y1 = flux[2] * (x/freq[2])**(-(p-1)/2)
p=2.5
y2 = flux[2] * (x/freq[2])**(-(p-1)/2)
plt.fill_between(x, y1, y2, color='grey')

#Plot nu**(-p/2)
x = np.linspace(freq[2], freq[-1])
p=2
y1 = flux[2] * (x/freq[2])**(-p/2)
p=2.5
y2 = flux[2] * (x/freq[2])**(-p/2)
plt.fill_between(x, y1, y2, color='grey')


plt.text(9, 170, "Spectrum at $\Delta t = 10\,$days", fontsize=14)

plt.xlabel("Frequency [GHz]", fontsize=16)
plt.ylabel("Flux [mJy]", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
#plt.show()

plt.savefig("spec_day_10.png")
