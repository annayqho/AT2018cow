""" Print table with VLASS limits """

import numpy as np
import sys
sys.path.append("/Users/annaho/Github/Query_VLASS")
from vlass_search import search_vlass
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/code")
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import ascii
from ret_radio import get_transients

def get_lims():
    headings = np.array(['ID', 'RA [hh:mm:ss]', 'Dec [dd:mm:ss]', '$\Delta t$', 'Limit ($\mu$Jy)'])
    label = "vlass"
    caption = "VLASS Limits for Rapidly Evolving Transients"
    names, ra_raw, dec_raw, dates, z = get_transients()
    dt = []
    limits = []

    for ii,name in enumerate(names):
        c = SkyCoord(ra_raw[ii], dec_raw[ii], unit=(u.hourangle, u.deg))
        out = search_vlass(name, c, dates[ii])
        if out is not None:
            lim, obsdate =out
            limits.append(lim)
            dt.append(obsdate-dates[ii])
        else:
            limits.append(None)
            dt.append(None)
    return names, ra_raw, dec_raw, dates, z, dt, limits


names, ra_raw, dec_raw, dates, z, dt, limits = get_lims()
ascii.write([names, ra_raw, dec_raw, dates, z, dt, limits], 'radio_lims.dat',
    names=['Name', 'RA', 'Dec', 'Date', 'z', 'dt', 'limit'])
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

for ii,ID in enumerate(names):
    if np.logical_or(limits[ii] == None, limits[ii] == 'nan'):
        limitstr = '-'
        tstr = '-'
    else:
        limitstr = int(limits[ii])
        tstr = int(dt[ii])
        row = rowstr %(ID, ra_raw[ii], dec_raw[ii], tstr, limitstr)
        # only print if not none
        outputf.write(row)
outputf.write("\enddata \n")
outputf.write("\end{deluxetable} \n")
outputf.close()
