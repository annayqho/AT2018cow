""" Search VLASS for a given RA and Dec """

import numpy as np
import os
from urllib.request import urlopen
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


def get_subtiles(tilename):
    """ For a given tile name, get the filenames in the VLASS directory.
    Parse those filenames and return a list of subtile RA and Dec.
    RA and Dec returned as a SkyCoord object
    """
    urlpath = urlopen('https://archive-new.nrao.edu/vlass/quicklook/VLASS1.1/%s/' %tilename)
    string = (urlpath.read().decode('utf-8')).split("\n")
    vals = np.array([val.strip() for val in string])
    keep_link = np.array(["href" in val.strip() for val in string])
    keep_name = np.array([tilename in val.strip() for val in string])
    string_keep = vals[np.logical_and(keep_link, keep_name)]
    fname = np.array([val.split("\"")[1] for val in string_keep])
    pos_raw = np.array([val.split(".")[4] for val in fname])
    ra_raw = np.array([val.split("-")[0] for val in pos_raw])
    dec_raw = np.array([val.split("-")[1] for val in pos_raw])
    ra = []
    dec = []
    for ii,val in enumerate(ra_raw):
        hms = "%sh%sm%ss" %(val[1:3], val[3:5], val[5:])
        ra.append(hms)
        dms = "%sd%sm%ss" %(
                dec_raw[ii][0:2], dec_raw[ii][2:4], dec_raw[ii][4:])
        dec.append(dms)
    ra = np.array(ra)
    dec = np.array(dec)
    c = SkyCoord(ra, dec, frame='icrs')
    return fname, c



if __name__=="__main__":
    # Take the lowest-redshift object
    ra = '02h45m06.67s'
    dec = '-00d44m42.77s'

    # Convert dec to degrees and RA to hours
    c = SkyCoord(ra, dec, frame='icrs')
    ra_h = c.ra.hour
    dec_d = c.dec.deg

    tiles = get_tiles()
    tilename = search_tiles(tiles, ra_h, dec_d)
    print(tilename)

    subtiles, c_tiles = get_subtiles(tilename)
    dist = c.separation(c_tiles)
    subtile = subtiles[np.argmin(dist)]

    url_get = "https://archive-new.nrao.edu/vlass/quicklook/VLASS1.1/%s/%s" %(
            tilename, subtile)
    fname = url_get + "%s.I.iter1.image.pbcor.tt0.subim.fits" %subtile[0:-1]
    print(fname)
    # and then wget this file


