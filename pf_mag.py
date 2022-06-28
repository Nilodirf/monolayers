import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as sp

from readinput import param

def itmag(ss, te, dmloc, mus, pf, pl):
    dmus=(param['rhosd']*dmloc-mus/pl)*param['dt']
    return(dmus)

def locmag(ss, mz, te, mus, fs, sup, sdn, pf, pl):
    const_ijk=1/pf*(param[`Jloc`]*mp*mz-param['J']*mz*mp)/np.sinh((param[`Jloc`]*mp*mz/param[`locspin`]-param['J']/param['s']*mz*mp)/(2*sp.k*param['tc']))
    fsup=sup*fs
    fsdn=sdn*fs
    wupwegprep=const_ijk0*np.exp(-(param['J']*mz-mus)/2/param['s']/sp.k/te[:ss[2]])
    wdnwegprep=const_ijk*np.exp((param['J']*mz-mus)/2/param['s']/sp.k/te[:ss[2]])
    wupweg=wupwegprep[...,np.newaxis]*fsup
    wdnweg=wdnwegprep[...,np.newaxis]*fsdn
    wuphin=np.roll(wupweg,1)
    wdnhin=np.roll(wdnweg,-1)
    dfs=(-wupweg-wdnweg+wuphin+wdnhin)*param['dt']
    return(dfs)