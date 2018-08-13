import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
import numpy as np

freq = np.array([15, 230, 353])
flux = np.array([0.5, 30.2, 33.02])
flux_err = np.array([0.1, 1.8, 0.46])

plt.errorbar(freq, flux, yerr=flux_err, fmt='.', c='k', marker='o')
plt.plot(freq, flux, linestyle='-', lw=0.5, c='k')

# Plot nu**2
x = np.linspace(freq[0], 102)
y = flux[0] * (x/freq[0])**2
plt.plot(x, y, linestyle='--', lw=2.0, c='k')
plt.text(30, 6, "$F_\\nu \propto \\nu^2$", fontsize=14)

# Plot nu**(1/3)
x = np.linspace(102, freq[1])
y = flux[1] * (x/freq[1])**(1/3)
plt.plot(x, y, linestyle='--', lw=2.0, c='k')
plt.text(113, 31, "$F_\\nu \propto \\nu^{1/3}$", fontsize=14)

#Plot nu**(-(p-1)/2)
# maybe not applicable in this case
#x = np.linspace(freq[2], 400)
#p=2
#y1 = flux[2] * (x/freq[2])**(-(p-1)/2)
#p=2.5
#y2 = flux[2] * (x/freq[2])**(-(p-1)/2)
#plt.fill_between(x, y1, y2)

plt.text(14, 30, "Spectrum at $\Delta t = 6\,$days", fontsize=14)

plt.xlabel("Frequency [GHz]", fontsize=16)
plt.ylabel("Flux [mJy]", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
#plt.show()

plt.savefig("spec_day_6.png")
