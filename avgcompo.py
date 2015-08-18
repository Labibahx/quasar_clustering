"""13 Aug 2015 
read names of fits files from fits files and create a mean composite spectrum from them
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.stats import sigmaclip