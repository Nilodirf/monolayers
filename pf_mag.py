import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as sp

from readinput import param

def itmag(ss, te, dmloc, magp, pf, pl):
    itmom = param['muat'] - param['locmom']
    dmagp=-param['locmom']/itmom*dmloc-magp/pl*param['dt']
    return(dmagp)

def locmag(ss, mz, te, mp, fs, sup, sdn, pf, pl):
    print(fs)
    itmom = param['muat'] - param['locmom']
    p_coup= 1.6e-19
    f_coup= 50e-19*abs(mz)
    const_ijk=1/pf*(f_coup-p_coup)/sp.k/param['tc']/np.sinh((f_coup-p_coup)/(2*sp.k*te[:ss[2]]))
    fsup=sup*fs
    fsdn=sdn*fs
    wupwegprep=const_ijk*np.exp(-(f_coup-p_coup)/(2*sp.k*te[:ss[2]]))*(mp+1)/2
    wdnwegprep=const_ijk*np.exp((f_coup-p_coup)/(2*sp.k*te[:ss[2]]))*(1-mp)/2
    wupweg=wupwegprep[...,np.newaxis]*fsup
    wdnweg=wdnwegprep[...,np.newaxis]*fsdn
    wuphin=np.roll(wupweg,1)
    wdnhin=np.roll(wdnweg,-1)
    dfs=(-wupweg-wdnweg+wuphin+wdnhin)*param['dt']
    return(dfs)