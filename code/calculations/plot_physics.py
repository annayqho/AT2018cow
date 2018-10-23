""" Measure the various physical quantities on each day """

import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import Planck15
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/code")
from get_radio import get_spectrum
from synchrotron_fit import self_abs, fit_self_abs


# CONSTANTS
c = 3E10
m_p = 1.67E-24
c6 = 8.16E-41
c5 = 6.29E-24
c1 = 6.27E18
El = 8.17E-7
f = 0.5
q_e = 4.8E-10
m_e = 9.1E-28
sigmat = 6.6524E-25
eps_e = 0.1
eps_B = 0.01


def get_R_full(Fp, nup):
    """ This is Equation 11 in Chevalier 98 """
    top = 6 * c6**(gamma+5) * Fp**(gamma+6) * d**(2*gamma+12)
    bottom = alpha*f*(gamma-2)*np.pi**(gamma+5)*c5**(gamma+6)*El**(gamma-2)
    Rp = (top/bottom)**(1/(2*gamma+13)) * (nup/(2*c1))**(-1)
    return Rp


def get_R(Fp, nup, d_mpc):
    """ This is Equation 13 in C98
    Fp in Jy
    nup in GHz

    Returns Rp in cm
    """
    alpha = eps_e / eps_B
    Rp = 8.8E15 * alpha**(-1/19) * (f/0.5)**(-1/19) * (Fp)**(9/19) * \
         d_mpc**(18/19) * (nup/5)**(-1)
    return Rp


def get_B(Fp, nup, d_mpc):
    """ This is Equation 14 in C98
    Fp in Jy
    nup in GHz
    
    Returns Bp in Gauss
    """
    alpha = eps_e / eps_B
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


def run(dt, nupeak, fpeak, d, p, d_mpc):
    """ nupeak in GHz, fpeak in Jy, dist in cm """
    print("lpeak", fpeak*1E-23*4*np.pi*d**2)
    R = get_R(fpeak, nupeak, d_mpc)
    print("R", R/10**15)
    B = get_B(fpeak, nupeak, d_mpc)
    print("B", B)
    V = (4/3) * f * np.pi * R**3
    v = R/(86400*dt)
    beta = v/c
    print("Beta", beta)
    P = B**2/(8*np.pi*eps_B)
    rho = (4*P/3)/v**2
    print("rho", rho)
    n_p = rho/m_p
    n_e = n_p
    print("ne", n_e)
    UB = (B**2)/(8*np.pi) * V
    print("E", UB/eps_B)
    tauff = 8.235E-2 * n_e**2 * (R/3.086E18) * (8000)**(-1.35)
    print("tau_ff", tauff)
    Te = 1.2E6
    nuff = (1/(68*(Te/8000)**(-1.35)))**(-1/2.1)
    print("nu_ff [GHz]", nuff)
    Lff = 1.43E-27 * n_e**2 * Te**(1/2) * (4/3) * np.pi * (6E16)**3
    print("L_ff", Lff)
    # Gamma_m
    gammam = eps_e * (m_p / m_e) * ((p-2)/(p-1)) * beta**2
    print("gamma_m", gammam)
    # Gyrofrequency
    gammag = q_e*B / (2 * np.pi * m_e * c)
    print("gyrofrequency [MHz]", gammag/1E6)
    # nu_m
    num = gammam**2 * gammag
    print("nu_m [GHz]", num/1E9)
    # Cooling frequency
    gammac = 6 * np.pi * m_e * c / (sigmat * B**2 * dt * 86400)
    print("gammac", gammac)
    nuc = gammac**2 * q_e * B / (2 * np.pi * m_e * c)
    print("nu_cc [GHz]", nuc/1E9)


def at2018cow():
    nupeak = 100
    fpeak = 94E-3
    p = 3
    d = Planck15.luminosity_distance(z=0.014).cgs.value
    d_mpc = Planck15.luminosity_distance(z=0.014).value
    run(22, nupeak, fpeak, d, p, d_mpc)


def tde():
    """ we use the values at 4 different days """
    dt = 5+3.04
    nupeak = 600
    fpeak = 80E-3
    dt = 10+3.04
    nupeak = 250
    fpeak = 40E-3
    dt = 15+3.04
    nupeak = 140
    fpeak = 30E-3
    d = Planck15.luminosity_distance(z=0.354).cgs.value
    d_mpc = Planck15.luminosity_distance(z=0.354).value
    gamma = 3
    # need to add 3.04 to all days
    run(dt, nupeak, fpeak, d, gamma, d_mpc)


def sn2003L():
    """ Soderberg 2005
    "At 30 days, the peak flux density is 3.2 mJy at 22.5 GHz
    The observed optically thin spectral index is -1.1,
    same as for AT2018cow
    distance is 92 Mpc
    """
    d = Planck15.luminosity_distance(z=0.0205).cgs.value
    d_mpc = Planck15.luminosity_distance(z=0.0205).value
    beta = -1.1
    gamma = -2*beta+1
    nupeak = 22.5
    fpeak = 3.2E-3
    run(30, nupeak, fpeak, d, gamma, d_mpc)


def sn2007bg():
    """ Salas et al. 2013 
    In Phase 1 of the explosion, the peak luminosity reached
    is 4.1E28 erg/s/Hz observed on Day 55.9 at 8.46 GHz
    d = 152 Mpc
    """
    nupeak = 8.46
    fpeak = 1490E-6
    redshift = 0.0335
    d = Planck15.luminosity_distance(z=redshift).cgs.value
    d_mpc = Planck15.luminosity_distance(z=redshift).value
    gamma = 3.2 # they don't give it in the paper though
    dt = 55.9
    run(dt, nupeak, fpeak, d, gamma, d_mpc)


def sn2003bg():
    """ Soderberg 2006 
    Peak flux density of 85 mJy at 22.5 GHz, spectral index -1.1
    epoch is 35 days
    distance is 19.6 Mpc
    """
    redshift = 0.00442
    d = Planck15.luminosity_distance(z=redshift).cgs.value
    d_mpc = Planck15.luminosity_distance(z=redshift).value
    beta = -1.1
    gamma = -2*beta+1
    nupeak = 22.5
    fpeak = 85E-3
    run(35, nupeak, fpeak, d, gamma, d_mpc)


def sn1998bw():
    """ Kulkarni et al. 1998
    on Day 10 infer a peak frequency of 10 GHz and associated
    peak flux 50 mJy
    distance = 38 Mpc
    redshift is 0.0083
    optically thin spectral index is -0.75 on Day 40
    they use alpha = 1
    which I think means beta = -1, so that's 3.
    so I can still use the same constants
    """
    redshift = 0.0083
    d = Planck15.luminosity_distance(z=redshift).cgs.value
    d_mpc = Planck15.luminosity_distance(z=redshift).value
    beta = -1
    gamma = -2*beta+1
    nupeak = 10
    fpeak = 50E-3
    run(10, nupeak, fpeak, d, gamma, d_mpc)


def sn2009bb():
    """ Soderberg et al. 2010
    from their first spectrum at dt = 20,
    they infer nupeak = 6 GHz and a spectral peak luminosity
    of 3.6E28 erg/s/Hz
    distance is 40 Mpc
    spectral index is around p=3 as well
    """
    redshift = 0.009
    d = Planck15.luminosity_distance(z=redshift).cgs.value
    d_mpc = Planck15.luminosity_distance(z=redshift).value
    gamma = 3
    nupeak = 6
    fpeak = 1E23 * 3.6E28 / (4 * np.pi * d**2) 
    run(20, nupeak, fpeak, d, gamma, d_mpc)


def sn2006aj():
    """ Soderberg et al. 2006
    distance is 145 Mpc
    Day 5: radio emission peaks between 1.4 GHz and 4.9 GHz
    In the caption they say that at 5d the radio spectrum peaks near 4 GHz
    They don't give a peak flux but the flux of 4.86 GHz at 5d is 328 uJy
    so let's use that
    They assume a power-law index of 2.1...
    but let's assume 3 to make it consistent with the rest of the sources
    """
    redshift = 0.032
    d = Planck15.luminosity_distance(z=redshift).cgs.value
    d_mpc = Planck15.luminosity_distance(z=redshift).value
    gamma = 3
    nupeak = 4
    fpeak = 328E-6
    run(5, nupeak, fpeak, d, gamma, d_mpc)


def sn2010bh():
    """ Margutti 2013
    distance is z=0.0593
    no idea what the spectral index is
    so let's use gamma = 3
    they find that at 30 days, self-absorption frequency is around 5 GHz
    and the flux is 130 uJy
    """
    redshift = 0.0593
    d = Planck15.luminosity_distance(z=redshift).cgs.value
    d_mpc = Planck15.luminosity_distance(z=redshift).value
    gamma = 3
    nupeak = 5
    fpeak = 130E-6
    run(30, nupeak, fpeak, d, gamma, d_mpc)


def ptf11qcj():
    """ Corsi 2014
    distance is z=0.0287
    gamma = 3
    the peak luminosity at 5 GHz is around 7E28 erg/s/Hz
    and that was at 10 days
    """
    redshift = 0.0287
    d = Planck15.luminosity_distance(z=redshift).cgs.value
    d_mpc = Planck15.luminosity_distance(z=redshift).value
    gamma = 3
    nupeak = 5
    fpeak = 1E23 * 7E28 / (4 * np.pi * d**2) # in Jy
    run(10, nupeak, fpeak, d, gamma, d_mpc)


def sn88Z():
    """ van Dyk et al. 1993
    the redshift was z=0.022
    the 6 cm maximum flux density was 1.90 mJy
    the spectral index is alpha = -0.7 
    they get 9E32 erg/s as the luminosity
    that was at 1253 days after the explosion
    but I think that means that p = 2.4, not 3.2
    TO DO: look up what the c constants are for that case
    """
    redshift = 0.022
    d = Planck15.luminosity_distance(z=redshift).cgs.value
    d_mpc = Planck15.luminosity_distance(z=redshift).value
    gamma = 2.4
    nupeak = 5
    fpeak = 1.90E-3
    run(1253, nupeak, fpeak, d, gamma, d_mpc)
    

def sn79C():
    """ Weiler 
    I'm going to use the peak of the 1.4 GHz light curve,
    which is what I have in the mm light curve plot.
    nu = 1.4E9 
    it's around 12 mJy at 1400 days
    """
    d = 5.341805643483106e+25
    d_mpc = d/3E24
    nupeak = 1.4
    fpeak = 12E-3
    gamma = 3 # no idea
    run(1400, nupeak, fpeak, d, gamma, d_mpc)


def grb031203():
    """
    Soderberg
    the peak luminosity is 1E29 erg/s/Hz at 8.5 GHz
    and the peak time of that is 35 days or so
    the derived synchrotron parameters at 1 day:
    3.2E8 Hz is nu_a
    F_nu,a = 0.04 mJy
    """
    z = 0.105
    d = Planck15.luminosity_distance(z=0.105).cgs.value
    d_mpc = Planck15.luminosity_distance(z=0.105).value
    nupeak = 0.328
    fpeak = 0.04E-3
    gamma = 3 # no idea
    run(1, nupeak, fpeak, d, gamma, d_mpc)
    

if __name__=="__main__":
    at2018cow()
