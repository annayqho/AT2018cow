""" Search VLASS for a given RA and Dec """

import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord


def get_tiles():
    """ Get tiles 
    I ran wget https://archive-new.nrao.edu/vlass/VLASS_dyn_summary.php
    """

    inputf = open("VLASS_dyn_summary.php", "r")
    lines = inputf.readlines()
    inputf.close()

    header = list(filter(None, lines[0].split("  ")))
    # get rid of white spaces
    header = np.array([val.strip() for val in header])

    names = []
    dec_min = []
    dec_max = []
    ra_min = []
    ra_max = []
    obsdate = []

    # Starting at lines[3], read in values
    for line in lines[3:]:
        dat = list(filter(None, line.split("  "))) 
        dat = np.array([val.strip() for val in dat]) 
        names.append(dat[0])
        dec_min.append(float(dat[1]))
        dec_max.append(float(dat[2]))
        ra_min.append(float(dat[3]))
        ra_max.append(float(dat[4]))
        obsdate.append(dat[6])

    names = np.array(names)
    dec_min = np.array(dec_min)
    dec_max = np.array(dec_max)
    ra_min = np.array(ra_min)
    ra_max = np.array(ra_max)
    obsdate = np.array(obsdate)

    return (names, dec_min, dec_max, ra_min, ra_max, obsdate)


def search_tiles(tiles, ra_h, dec_d):
    """ Now that you've processed the file, search for the given RA and Dec """
    names, dec_min, dec_max, ra_min, ra_max, obsdate = tiles
    has_dec = np.logical_and(dec_d > dec_min, dec_d < dec_max)
    has_ra = np.logical_and(ra_h > ra_min, ra_h < ra_max)
    in_tile = np.logical_and(has_ra, has_dec)
    name = names[in_tile]
    if len(name) > 1:
        print("Error: this source is in more than one tile")
    else:
        return name[0]


if __name__=="__main__":
    # Take the lowest-redshift object
    ra = '02h45m06.67s'
    dec = '-00d44m42.77s'

    # Convert dec to degrees and RA to hours
    c = SkyCoord(ra, dec, frame='icrs')
    ra_h = c.ra.hour
    dec_d = c.dec.deg

    tiles = get_tiles()
    name = search_tiles(tiles, ra_h, dec_d)
    print(name)
