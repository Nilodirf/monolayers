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
for i in par['sdrate']:
    for j in par['sdissrate']:
        output.output(par,i, j)


