""" Print table with VLASS limits """

import numpy as np
from vlass_search import search_vlass
from ret_radio import get_transients

headings = np.array(['ID', 'RA', 'Dec', '$\Delta t$ (d)', 'Limit ($\\mu$Jy)'])
label = "vlass"
caption = "VLASS Limits for Rapidly Evolving Transients"
names, ra_raw, dec_raw, dates = get_transients()
limits, obsdates = search_vlass()

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
    if np.logical_or(dates[ii] < 0, limits[ii] == 'nan'):
        datestr = '-'
        limitstr = '-'
    else:
        datediff = obsdates[ii] - dates[ii]
        if datediff < 0:
            datestr = '-'
            limitstr = '-'
        else:
            datestr = str(datediff)
            limitstr = str(limits[ii])
    row = rowstr %(ID, ra_raw[ii], dec_raw[ii], datestr, limitstr)
    print(row)
    outputf.write(row)
outputf.write("\enddata \n")
outputf.write("\end{deluxetable} \n")
outputf.close()
