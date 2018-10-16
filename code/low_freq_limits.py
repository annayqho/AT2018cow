""" Get low-frequency VLASS limits """

from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from get_radio import get_data_all


# Step 1: get all of the days with ATCA data
tel, freq, days, flux, eflux_form, eflux_sys = get_data_all()
choose = np.logical_and(tel=='ATCA')

# Step 2: for each of those days, fit nu squared and extrapolate
# down to 3 GHz, to generate a 3 GHz light curve

# Step 3: Plot the 3 GHz light curve as luminosity
