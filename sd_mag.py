import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as sp

from readinput import param

def sd_mag(ts, ss, mz, te, mus):
    dmz=1/param['tausd']*(mz-(mus/(2*sp.k*param['tc'])))*(1-mz/np.tanh(2*mz*sp.k*param['tc']-mus/(2*sp.k*te[:ss[2]])))*param['dt']
    dmus=param['rhosd']*dmz-mus/param['taus']
    return(dmz, dmus)