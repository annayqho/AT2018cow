""" Make a grid of light curves from the SMA, ATCA, and Swift/XRT """

import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
import numpy as np

days = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

# flux at 231.5 GHz
