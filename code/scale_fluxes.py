""" Scale fluxes using variations in 1635 at that frequency

Notes that I deleted Jul 23 from the file and June 28
and deleted comments from July 2
"""

import numpy as np
from astropy.time import Time

def sma_lc():
    data_dir = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
    lines = open("%s/SMA_AT2018cow_quasar.txt" %data_dir, "r").readlines()
    t0 = Time('2018-06-16') # MJD 58285.0, what Dan has in his light curves
    dt = []
    nu = []
    f = []
    ef = []
    ref = []

    for line in lines:
        has_date = 'Jun' in line or 'Jul' in line or 'Aug' in line
        if has_date and len(line.split()) > 10:
            date = line.split()[1]
            if date == 'Jun':
                mm = '06'
            elif date == 'Jul':
                mm = '07'
            elif date == 'Aug':
                mm = '08'
            else:
                print("something went wrong")

            dt.append((Time('2018-%s-%s' %(mm,line.split()[2]))-t0).value)
            nu.append(line.split()[6])
            f.append(line.split()[7])
            ef.append(line.split()[9])
            ref.append(line.split()[10])
    dt = np.array(dt, dtype=float)
    nu = np.array(nu, dtype=float)
    f = np.array(f, dtype=float)
    ef = np.array(ef, dtype=float)
    ref = np.array(ref, dtype=float)

    # Step 1: Measure the average quasar flux at each frequency
    u_nu = np.unique(nu)
    u_ref = np.array([np.mean(ref[nu==val]) for val in u_nu])

    # Step 2: Scale each flux measurement by the 
    # ratio of this avg/1635 flux at that epoch
    npts = len(dt)
    f_scaled = np.zeros(npts)
    ef_scaled = np.zeros(npts)

    for ii,fmeas in enumerate(f):
        mean_ref = u_ref[u_nu==nu[ii]][0]
        ratio = mean_ref / ref[ii]
        f_scaled[ii] = f[ii] * ratio
        ef_scaled[ii] = ef[ii] * ratio

    # Generate the sma_radio_lc_231.dat file
    # This is the light curve of 231.5 GHz
    dt_final = []
    flux_final = []
    eflux_final = []

    for ii,dt_val in enumerate(np.unique(dt)):
        dt_final.append(dt_val)
        # Count the number of points on this day measured at 231.5 GHz
        choose = np.logical_and(nu[dt==dt_val] <= 234.6, nu[dt==dt_val] >= 230.6)
        # If the answer is 1, then store that point
        if sum(choose) == 1:
            flux_final.append(f_scaled[dt==dt_val][choose][0])
            eflux_final.append(ef_scaled[dt==dt_val][choose][0])
        # If the answer is 2, then take the average of those points
        if sum(choose) == 2:
            w = 1/(ef_scaled[dt==dt_val][choose])**2
            fmean, wsum = np.average(
                    f_scaled[dt==dt_val][choose], 
                    weights=w,
                    returned=True)
            efmean = np.sqrt(1/np.sum(w))
            flux_final.append(fmean)
            eflux_final.append(efmean)
        # If the answer is > 2, something is wrong...
        if sum(choose) > 2:
            print("something is wrong!")
        if sum(choose) == 0:
            # the only days that aren't within this tiny range
            # have nu = 231.5 and nu = 218.0,
            # so these are the ones I will bother to scale
            # choose the points with 243.3 and scale them to 231.5
            toscale = nu[dt==dt_val] == 243.3
            if sum(toscale) == 1:
                flux_final.append(f_scaled[dt==dt_val][toscale][0] * (243.3/231.5))
                eflux_final.append(
                        ef_scaled[dt==dt_val][toscale][0] * (243.3/231.5))
            else:
                toscale = nu[dt==dt_val] == 218
                if sum(toscale) == 1:
                    flux_final.append(f_scaled[dt==dt_val][toscale][0] * (218.0/231.5))
                    eflux_final.append(
                            ef_scaled[dt==dt_val][toscale][0] * (218.0/231.5))

    return np.array(dt_final), np.array(flux_final), np.array(eflux_final)
