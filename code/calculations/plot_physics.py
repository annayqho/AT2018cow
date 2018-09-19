""" Measure the various physical quantities on each day """

import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import Planck15
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/code")
from get_radio import get_spectrum
from synchrotron_fit import self_abs, fit_self_abs


d = Planck15.luminosity_distance(z=0.014).cgs.value
d_mpc = Planck15.luminosity_distance(z=0.014).value
c = 3E10
mp = 1.67E-24
c6 = 8.16E-41
c5 = 6.29E-24
c1 = 6.27E18
El = 8.17E-7
gamma = 3.4
alpha = 1
f = 0.5


def get_R_full(Fp, nup):
    """ This is Equation 11 in Chevalier 98 """
    top = 6 * c6**(gamma+5) * Fp**(gamma+6) * d**(2*gamma+12)
    bottom = alpha*f*(gamma-2)*np.pi**(gamma+5)*c5**(gamma+6)*El**(gamma-2)
    Rp = (top/bottom)**(1/(2*gamma+13)) * (nup/(2*c1))**(-1)
    return Rp


def get_R(Fp, nup):
    """ This is Equation 13 in C98
    Fp in Jy
    nup in GHz

    Returns Rp in cm
    """
    Rp = 8.8E15 * alpha**(-1/19) * (f/0.5)**(-1/19) * (Fp)**(9/19) * \
         d_mpc**(18/19) * (nup/5)**(-1)
    return Rp


def get_B(Fp, nup):
    """ This is Equation 14 in C98
    Fp in Jy
    nup in GHz
    
    Returns Bp in Gauss
    """
    Bp = 0.58 * alpha**(-4/19) * (f/0.5)**(-4/19) * (Fp)**(-2/19) * \
         d_mpc**(-4/19) * (nup/5)
    return Bp
    


def multiple_days():
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


if __name__=="__main__":
    R = get_R(94E-3, 100)
    print(R/10**15)
    B = get_B(94E-3, 100)
    dt = 22
    V = (4/3) * np.pi * R**3
    v = R/(86400*dt)
    beta = v/c
    print(beta)
    UB = (B**2)/(8*np.pi) * V

