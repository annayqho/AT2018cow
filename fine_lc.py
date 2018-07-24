""" Light curves binned into smaller time intervals """

import numpy as np
import matplotlib.pyplot as plt

def weighted_mean(flux, flux_err):
    mean = np.sum(flux/flux_err**2, axis=0) / np.sum(1/flux_err**2, axis=0)
    err_mean = np.sqrt(1 / np.sum(1/flux_err**2, axis=0))
    return mean, err_mean


# June 21, Day 5
raw = ["RxA LSB  11.31 +/-1.94  10.48 +/-1.13  15.44 +/-1.11  15.28 +/-1.27",
"RxA USB  13.07 +/-2.33  12.41 +/-1.33  14.96 +/-1.31  13.43 +/-1.53",
"RxB LSB  7.34 +/-2.39  12.72 +/-1.57  14.05 +/-1.35  15.36 +/-1.52",
"RxB USB  11.81 +/-2.63  15.01 +/-1.70  14.93 +/-1.45  15.95 +/-1.64"]

bands = []
flux = []
flux_err = []

for line in raw:
    bands.append(line.split("  ")[0])
    vals = line.split("  ")[1:]
    flux.append(np.array(
            [float(val.split("+/-")[0]) for val in vals]))
    flux_err.append(np.array(
            [float(val.split("+/-")[1]) for val in vals]))
   
bands = np.array(bands)
sb = np.array([val.split(" ")[1] for val in bands])
flux = np.array(flux)
flux_err = np.array(flux_err)

choose = sb == 'LSB'
mean_flux, emean_flux = weighted_mean(flux[choose], flux_err[choose])
