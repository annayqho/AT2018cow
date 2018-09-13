""" Measure the various physical quantities on each day """

import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import Planck15
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/code")
from get_radio import get_spectrum
from synchrotron_fit import self_abs, fit_self_abs


d = Planck15.luminosity_distance(z=0.014).cgs.value
c = 3E10
mp = 1.67E-24

for dt in np.arange(8,31):
    nu,f = get_spectrum(dt)
    plt.scatter(nu, f, c='k')
    # For points < 40 GHz, fit a nu**2 power law
    choose = nu < 40
    order = np.argsort(nu[choose])
    ba = fit_self_abs(nu[choose][order],f[choose][order])
    ma = 2
    xfit = np.linspace(5, 200)
    yfit = 10**(ma*np.log10(xfit)+ba)
    plt.plot(xfit, yfit, ls=':', c='k')

    # For points > 200 GHz, fit whatever power law
    choose = nu > 200
    m,b = np.polyfit(np.log10(nu[choose]), np.log10(f[choose]), deg=1)
    xfit = np.linspace(50, 1000)
    yfit = 10**(m*np.log10(xfit)+b)
    plt.plot(xfit, yfit, ls=':', c='k')

    # Find the intersection of these two lines
    xin = (ba-b) / (m-ma)
    yin = m*xin+b

    # Calculations
    nu_a = 10**xin
    F_a = 10**yin
    alpha = m
    theta = 120 * 60**(-1/17) * F_a**(8/17) * nu_a**((2+alpha-35)/34)
    R = theta * 10**(-6) * (1/206265) * d
    V = (4/3) * np.pi * R**3
    v = R/(86400*dt)
    beta = v/c
    B = 20.8*(nu_a*F_a / ((R/10**15)**3))**(2/7)
    U = 2 * (B**2)/(8*np.pi) * V

    #8.4E45 * B**2 * (R/10**15)**3
    n = (B**2 / (8 * np.pi * mp * v**2))

    # Plot results
    R_plt = r"%s \times 10^{15}" %np.round(R/10**15,1)
    U_plt = r"%s \times 10^{50}" %np.round(U/10**50,1)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Frequency [GHz]")
    plt.ylabel("Flux [mJy]")
    plt.text(5,100,
    r"$R=%s$ cm, $\beta=%s$" %(R_plt, np.round(beta,2)))
    plt.text(5,70,
    r"$B=%s$ G, $U=%s$" %(int(B), U_plt))
    plt.text(5, 170, "Day %s" %dt, fontsize=15)
    plt.text(200, 120,
            r"$\nu_a=%s$ GHz" %(int(nu_a)))
    plt.text(200, 90,
            r"$F_a=%s$mJy" %(int(F_a)))
    plt.savefig("day_%s.png" %dt)
    plt.close()
