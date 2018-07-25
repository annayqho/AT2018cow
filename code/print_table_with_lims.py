""" Print table with VLASS limits """

import numpy as np
from vlass_search import search_vlass
from ret_radio import get_transients

headings = np.array(['ID', 'RA', 'Dec', 'Expl Date', 'Limit (uJy)'])
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
    if dates[ii] < 0:
        dates[ii] = '-'
    if limits[ii] == 'nan':
        limits[ii] = '-'
    row = rowstr %(ID, ra_raw[ii], dec_raw[ii], dates[ii], limits[ii])
    outputf.write(row)
outputf.write("\enddata \n")
outputf.write("\end{deluxetable} \n")
outputf.close()
