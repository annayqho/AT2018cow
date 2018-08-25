""" Print table with VLASS limits """

import numpy as np
import sys
sys.path.append("/Users/annaho/Github/Query_VLASS")
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/code")
from astropy.coordinates import SkyCoord
from astropy import units as u
from vlass_search import search_vlass
from ret_radio import get_transients

headings = np.array(['ID', 'RA', 'Dec', 'Expl Date', 'Limit'])
label = "vlass"
caption = "VLASS Limits for Rapidly Evolving Transients"
names, ra_raw, dec_raw, dates, z = get_transients()
limits = []

for ii,name in enumerate(names):
    c = SkyCoord(ra_raw[ii], dec_raw[ii], unit=(u.hourangle, u.deg))
    limits.append(search_vlass(name, c, dates[ii]))

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
        limitstr = limits[ii]
        tstr = dates[ii]
    row = rowstr %(ID, ra_raw[ii], dec_raw[ii], tstr, limitstr)
    outputf.write(row)
outputf.write("\enddata \n")
outputf.write("\end{deluxetable} \n")
outputf.close()
