import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as sp

from readinput import param

def itmag(ss, te, dmloc, mus):
    dmus=param['rhosd']*dmloc-mus/param['sdissrate']
    return(dmus)

def locmag(ss, mz, te, mus, fs, sup, sdn):
    const_ijk=param['sdrate']*(param['J']*mz/param['s']-mus)/np.sinh((param['J']*mz-mus)/2/param['s']/sp.k/te[:ss[2]])
    fsup=sup*fs
    fsdn=sdn*fs
    wupwegprep=const_ijk*np.exp(-(param['J']*mz-mus)/2/param['s']/sp.k/te[:ss[2]])
    wdnwegprep=const_ijk*np.exp((param['J']*mz-mus)/2/param['s']/sp.k/te[:ss[2]])
    wupweg=wupwegprep[...,np.newaxis]*fsup
    wdnweg=wdnwegprep[...,np.newaxis]*fsdn
    wuphin=np.roll(wupweg,1)
    wdnhin=np.roll(wdnweg,-1)
    dfs=(-wupweg-wdnweg+wuphin+wdnhin)*param['dt']
    return(dfs)