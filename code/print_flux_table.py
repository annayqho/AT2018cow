""" Print table with VLASS limits """

import numpy as np
import sys
sys.path.append("/Users/annaho/Github/Query_VLASS")
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/code")
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from scale_fluxes import sma_lc

headings = np.array(
        ['$\Delta t_\mathrm{expl}$ (d)', 'Facility', 
         'Frequency (GHz)', 'Flux Density (mJy)'])
label = "flux"
caption = "Flux measurements for AT2018cow. \
        SMA measurements have formal uncertainties shown, \
        which are appropriate for in-band measurements on a given night. \
        However, for night-to-night comparisons, true errors are dominated by \
        systematics and are roughly 10\%--15\% unless indicated otherwise. \
        ALMA measurements have roughly 5\% uncertainties in Bands 3 and 4, \
        10\% uncertainties in Bands 7 and 8, and a 20\% uncertainty in Band 9.\
        ATCA measurements have formal errors listed, \
        but also have systematic uncertainties of roughly 10\%."

# Print the table headers
ncol = len(headings) 
colstr = ""
colstr += 'l'
for col in np.arange(ncol-1): colstr+="r"
print(colstr)

colheadstr = ""
for col in np.arange(ncol-1):
    colheadstr += "\colhead{%s} & " %headings[col]
colheadstr += "\colhead{%s}" %headings[-1]

rowstr = ""
for col in np.arange(ncol-1):
    rowstr += "%s & " 
rowstr += "%s \\\ \n"

outputf = open("paper_table_%s.txt" %label, "w")
outputf.write("\\startlongtable \n")
outputf.write("\\begin{deluxetable}{%s} \n" %colstr)
outputf.write("\\tablecaption{%s\label{tab:%s}} \n" %(caption,label))
outputf.write("\\tablewidth{0pt} \n")
outputf.write("\\tablehead{ %s } \n" %colheadstr)
#outputf.write("\\rotate \n")
outputf.write("\\tabletypesize{\scriptsize} \n")
outputf.write("\startdata \n")


# Get the SMA light curves with fractional dt and scaled flux
a,b,c = sma_lc()
dt, nu, f, ef = a
tvals = np.unique(dt)

for ii,t in enumerate(tvals):
    # For each time
    choose = dt == t
    # For each frequency observed at that time 
    nu_unique = np.unique(nu[choose])
    for jj,nuval in enumerate(nu_unique):
        # Select the relevant data points
        flux_raw = f[choose][nu[choose] == nuval]
        eflux_raw = ef[choose][nu[choose] == nuval]
        if len(flux_raw) == 2:
            w = 1/(eflux_raw)**2
            fmean, wsum = np.average(flux_raw, weights=w, returned=True)
            efmean = np.sqrt(1/np.sum(w))
            fstr = '$%s \pm %s$' %(np.round(fmean,2), np.round(efmean,2))
        elif len(flux_raw) == 1:
            fstr = '$%s \pm %s$' %(np.round(flux_raw[0],2), np.round(eflux_raw[0],2))
        else:
            print("something is wrong")
        if int(t) == 31:
            fstr += '$^a$'
        row = rowstr %(np.round(t,2), 'SMA', nuval, fstr)
        outputf.write(row)

# Next, get the ATCA light curves with fractional dt and scaled flux
# Day 10: 09:00 - 14:00 UT, so mean=11.30, which is fractional day 0.48
# Day 13: 08:30 - 14:00 UT, so mean=11:15, which is a fractional day 0.47
# Day 17: 09:00 - 13:45 UT, so mean=11:22.5, which is fractional day 0.47
dt = np.array(
        [10.48, 10.48, 10.48, 13.47, 13.47, 13.47, 13.47, 13.47, 17.47, 17.47, 
         19.615, 28.44, 34.43, 81.37])
nu = np.array(
        [5.5, 9, 34, 5.5, 9, 16.7, 21.2, 34, 5.5, 9, 
         34, 34, 34, 34])
f = np.array(
        [0.15, 0.27, 5.6, 0.22, 0.52, 1.5, 2.3, 7.6, 0.41, 0.99, 
         14.26, 30.59, 42.68, 6.97])
ef = np.array(
        [99, 0.06, 0.16, 0.05, 0.04, 0.1, 0.3, 0.5, 0.04, 0.03, 
         0.21, 0.20, 0.19, 0.09])
for ii,t in enumerate(dt):
    if ii == 0:
        if nu[ii] == 5.5:
            fstr = '< %s' %f[ii]
        else:
            fstr = '$%s \pm %s$' %(f[ii], ef[ii])
    else:
        fstr = '$%s \pm %s$' %(f[ii], ef[ii])
    row = rowstr %(t, 'ATCA', nu[ii], fstr)
    outputf.write(row)


# Finally, get the ALMA light curves with fractional dt and scaled flux
# I went to this site: http://almascience.nrao.edu/aq/
# I assume that the times are in UT...?
# Yes, Bjorn confirms that they are in UT
t0 = Time(58285, format='mjd')
t = []
t.extend(['2018-06-30T00:48:10']*4) # Band 7
t.extend(['2018-06-30T03:22:38']*4) # Band 8
t.extend(['2018-07-08T00:35:42']*4) # Band 4
t.extend(['2018-07-08T00:53:03']*4) # Band 5
t.extend(['2018-07-09T01:31:41']) # Band 9
t = Time(t, format='isot', scale='utc')
dt = (t-t0).value
nu = np.array(
    [336.5, 338.5, 348.5, 350.5, 398, 400, 410, 412,
     90.5, 92.5, 102.5, 104.5, 138, 140, 150, 152,
     671])
f = np.array(
    [29.40, 29.10, 28.49, 28.29, 26.46, 26.21, 25.69, 25.95,
     91.18, 92.31, 93.97, 93.57, 85.10, 84.58, 80.62, 79.71,
     31.5])
ef = np.array(
    [2.94, 2.91, 2.85, 2.83, 2.65, 2.62, 2.57, 2.60, 
     0.46, 0.46, 0.47, 0.47, 0.43, 0.42, 0.40, 0.40,
     6.3])
for ii,t in enumerate(dt):
    fstr = '$%s \pm %s$' %(f[ii], ef[ii])
    row = rowstr %(np.round(t,2), 'ALMA', nu[ii], fstr)
    outputf.write(row)


outputf.write("\enddata \n")
outputf.write("\\tablecomments{\${^a}$ Systematic uncertainty 20\% due to uncertain flux calibration}")
outputf.write("\end{deluxetable} \n")
outputf.close()
