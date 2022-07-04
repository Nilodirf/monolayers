import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as sp

duration = 1000  # milliseconds
freq = 440  # Hz

#import other files
import readinput
import output

#read out inputfile
par=readinput.readout()
for power in par['pp']:
    output.output(par, power)


