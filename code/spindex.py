""" Analysis of the spectral index of AT2018cow """

import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt
import emcee
import corner
from get_radio import get_data_all, get_spectrum


def lnlike(theta, x, y, yerr):
    """ Log likelihood function """
    a, b, lnf = theta
    model = a*x**(b)
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))


def lnprior(theta):
    """ log-prior
    I think it's reasonable to have a constant prior
    across the distribution
    Allow the power to be anywhere between 0.1 and 5
    Allow the coefficient to be anywhere between 1 and 10000
    """
    a, b, lnf = theta
    if 5E2 < 5E7 and -3.0 < b < -0.1 and -10 < lnf < 1:
        return 0.0
    return -np.inf


def lnprob(theta, x, y, yerr):
    """ log-probability function """
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)


tel, freq, day, flux, eflux_form, eflux_sys = get_data_all()

choose_tel = np.logical_or(tel=='ALMA', tel=='SMA')
#choose_tel = tel == 'SMA'
choose = np.logical_and(choose_tel, day == 14)

nu = freq[choose]
print(nu)
f = flux[choose]
print(f)
ef = eflux_form[choose]
#ef[ef<0] = f[ef<0] * 0.0001
ef = flux[choose] * 0.001
#ef = eflux_sys[choose]
#ef = np.sqrt(eflux_form[choose]**2 + eflux_sys[choose]**2)
print(ef)

npts = len(nu)

# if there are fewer than 2 values, don't bother
if npts < 2:
    print("fewer than two points, can't calculate an index")

# if there are two values
if npts == 2:
    alpha = np.log10(f[-1]/f[0])/np.log10(nu[-1]/nu[0])
    ealpha = (1/np.log10(nu[-1]/nu[0])) * (ef[-1]/f[-1]-ef[0]/f[0])
    print(alpha, ealpha)
    xfit = np.linspace(min(nu), max(nu), 100)
    yfit = f[0] * 10**(alpha*np.log10(xfit/nu[0]))
    plt.scatter(nu, f, c='k')
    plt.plot(xfit, yfit, ls=':', c='r')
    plt.show()

# if there are more than two values, do an MCMC fit to the data
if npts > 2:
    # initialize a guess
    b_true = np.log10(f[3]/f[0])/np.log10(nu[3]/nu[0])
    a_true = f[1]/(nu[1]**b_true)
    print(a_true, b_true)
    f_true = 1.5

    ndim, nwalkers = 3, 100
    nll = lambda *args: -lnlike(*args)
    result = op.minimize(
            nll, [a_true, b_true, np.log(f_true)], args=(nu, f, ef))
    m_ml, b_ml, lnf_ml = result["x"]
    pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(nu, f, ef))
    sampler.run_mcmc(pos, 500)
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    #fig = corner.corner(samples, labels=["$a$", "$b$", "$\ln\,f$"])

    xl = np.linspace(200, 500, 1000)
    for a, b, lnf in samples[np.random.randint(len(samples), size=100)]:
        plt.plot(xl, a*xl**b, color="k", alpha=0.1)
    #plt.plot(xl, a_true*xl**b_true, color="r", lw=2, alpha=0.8)
    plt.errorbar(nu, f, yerr=ef, fmt=".k")

    samples[:, 2] = np.exp(samples[:, 2])
    a_mcmc, b_mcmc, f_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                 zip(*np.percentile(samples, [16, 50, 84],
                                                    axis=0)))
    print(b_mcmc)

#plt.show()
plt.savefig("mcmc_fit_day_14.png")
