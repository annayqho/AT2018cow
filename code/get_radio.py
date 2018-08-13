""" Get the radio data """


def get_data():
    """ Get the full set of data: 

    Returns
    -------
    tel: array telling you which telescope the data point is from
    freq: freq of the data point
    days: day of the data point, from time of explosion
    flux: flux of the data point
    flux_err: uncertainty in flux of the data point
    """
    data_dir = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
    dat = Table.read(
        "%s/radio_lc.dat" %data_dir, delimiter="&", format='ascii.no_header')
    days = np.array(dat['col1'])
    tel = np.array(dat['col2'])
    freq = np.array(dat['col3']).astype(float)
    flux_raw = np.array(dat['col4'])
    flux = np.zeros(len(flux_raw))
    flux_err = np.zeros(len(flux_raw))

    # Now, go through and apply the appropriate uncertainties
    for ii,val in enumerate(flux_raw):
        if tel[ii] == 'ALMA':
            flux[ii] = val.split("$")[1]
            if freq[ii] < 200: # bands 3 and 4
                flux_err[ii] = 0.05*flux[ii]
            elif freq[ii] > 600: # band 9
                flux_err[ii] = 0.20*flux[ii]
            else: # bands 7 and 8
                flux_err[ii] = 0.10*flux[ii]
        elif tel[ii] == 'SMA':
            flux[ii] = float(val.split("pm")[0][1:])
            flux_err_temp = float(val.split("pm")[1][0:-1])
            # add 10% uncertainty
            flux_err[ii] = np.sqrt(
                    (0.10*flux[ii])**2 + flux_err_temp**2)
        elif tel[ii] == 'ATCA':
            flux[ii] = float(val.split("pm")[0][1:])
            # add 10% uncertainty
            flux_err_temp = float(val.split("pm")[1][0:-1])
            if flux_err_temp < 0:
                flux_err[ii] = -99
            else:
                flux_err[ii] = np.sqrt(
                        (0.10*flux[ii])**2 + flux_err_temp**2)
    return tel, freq, days, flux, flux_err
