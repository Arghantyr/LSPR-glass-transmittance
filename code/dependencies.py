# Import dependencies

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from scipy import optimize
from scipy import integrate
from scipy import stats
import scipy.constants as constants

# Get the important physical constants
h = constants.physical_constants["Planck constant"][0]
c = constants.physical_constants["speed of light in vacuum"][0]
Jev = constants.physical_constants["joule-electron volt relationship"][0]
pi = np.pi
