""" Read the table from the radio data """

import numpy as np


DATA_DIREC = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data/radio_compilations"


def zauderer():
    lines = open(
            "%s/Zauderer2011/table.txt" %DATA_DIREC, "r").readlines()
    dt = []
    nu = []
    f = []
    islim = []
    ef = []

    for line in lines:
        dt.append(float(line.split("&")[1]))
        nu.append(float(line.split("&")[3]))
        if '<' in line:
            islim.append(True)
            f.append(float(
                line.split("&")[4].split("$")[1][1:]))
        else:
            islim.append(False)
            f.append(float(
                line.split("&")[4].split("$")[1].split("pm")[0]))
            ef.append(float(
                line.split("&")[4].split("$")[1].split("pm")[1]))

    dt = np.array(dt)
    nu = np.array(nu)
    f = np.array(f)
    ef = np.array(ef)
    islim = np.array(islim)

    return nu, dt, f, ef, islim


def read_2003L():
    """ this paper doesn't give upper limits """
    lines = open(
            "%s/SN2003L/table.txt" %DATA_DIREC, "r").readlines()
    dt = []
    nu = []
    f = []
    ef = []

    for line in lines:
        if 'nodata' not in str(line.split("&")[2]):
            dt.append(float(line.split("&")[1]))
            nu.append(4.9E9)
            f.append(float(
                line.split("&")[2].split("$")[0]))
            ef.append(float(
                line.split("&")[2].split("$")[2]))
        if 'nodata' not in str(line.split("&")[3]):
            dt.append(float(line.split("&")[1]))
            nu.append(8.5E9)
            f.append(float(
                line.split("&")[3].split("$")[0]))
            ef.append(float(
                line.split("&")[3].split("$")[2]))
        if 'nodata' not in str(line.split("&")[4]):
            dt.append(float(line.split("&")[1]))
            nu.append(15E9)
            f.append(float(
                line.split("&")[4].split("$")[0]))
            ef.append(float(
                line.split("&")[4].split("$")[2]))
        if 'nodata' not in str(line.split("&")[5]):
            dt.append(float(line.split("&")[1]))
            nu.append(22.5E9)
            f.append(float(
                line.split("&")[5].split("$")[0]))
            ef.append(float(
                line.split("&")[5].split("$")[2]))

    dt = np.array(dt)
    nu = np.array(nu)
    f = np.array(f)
    ef = np.array(ef)

    return nu, dt, f, ef


def line_1993J(dt, nu, f, ef, islim, line, col, wl):
    temp = str(line.split("&")[col])
    if 'nodata' not in temp:
        dt.append(float(line.split("&")[1]))
        nu.append(3E10/wl) # freq
        if '<' in temp:
            islim.append(True)
            f.append(float(
                temp.split("$")[1].split('<')[1]))
            ef.append(-99)
        else:
            islim.append(False)
            f.append(float(
                temp.split("$")[1].split("pm")[0]))
            ef.append(float(
                temp.split("$")[1].split("pm")[1]))
    return dt, nu, f, ef, islim


def read_1993J_low_freq():
    lines = open(
            "%s/SN1993J/table_low_freq.txt" %DATA_DIREC, "r").readlines()
    dt = []
    nu = []
    f = []
    ef = []
    islim = []

    for line in lines:
        dt, nu, f, ef, islim = line_1993J(dt, nu, f, ef, islim, line, 3, 20)
        dt, nu, f, ef, islim = line_1993J(dt, nu, f, ef, islim, line, 4, 6)
        dt, nu, f, ef, islim = line_1993J(dt, nu, f, ef, islim, line, 5, 3.6)
        dt, nu, f, ef, islim = line_1993J(dt, nu, f, ef, islim, line, 6, 2)
        dt, nu, f, ef, islim = line_1993J(dt, nu, f, ef, islim, line, 7, 1.2)

    nu = np.array(nu)
    f = np.array(f)
    ef = np.array(ef)
    dt = np.array(dt)
    islim = np.array(islim)

    return nu, dt, f, ef, islim


def read_1993J_high_freq():
    lines = open(
            "%s/SN1993J/table_high_freq.txt" %DATA_DIREC, "r").readlines()
    dt = []
    nu = []
    f = []
    ef = []
    islim = []

    for line in lines:
        dt.append(float(line.split("&")[1]))
        nu.append(1E9*float(line.split("&")[4]))
        if '<' not in line:
            islim.append(False)
            f.append(float(line.split("&")[3].split("$")[1].split("pm")[0]))
            ef.append(float(line.split("&")[3].split("$")[1].split("pm")[1]))
        else:
            islim.append(True)
            f.append(float(line.split("&")[3].split("$")[1].split("<")[1]))
            ef.append(-99)

    nu = np.array(nu)
    f = np.array(f)
    ef = np.array(ef)
    dt = np.array(dt)
    islim = np.array(islim)

    return nu, dt, f, ef, islim
