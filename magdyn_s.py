import numpy as np
from scipy import constants as sp

from readinput import param

def magdyn_s(ts, ss, mz, fs, sup, sdn, te, tp):
    
    const_ijk=param['asc']*tp[:ss[2]]*mz/4/param['s']/np.sinh(param['J']*mz/2/param['s']/sp.k/te[:ss[2]])*param['gepf'](te[:ss[2]])
    fsup=sup*fs
    fsdn=sdn*fs
    wupwegprep=const_ijk*np.exp(-param['J']*mz/2/param['s']/sp.k/te[:ss[2]])
    wdnwegprep=const_ijk*np.exp(param['J']*mz/2/param['s']/sp.k/te[:ss[2]])
    wupweg=wupwegprep[...,np.newaxis]*fsup
    wdnweg=wdnwegprep[...,np.newaxis]*fsdn
    wuphin=np.roll(wupweg,1)
    wdnhin=np.roll(wdnweg,-1)
    dfs=(-wupweg-wdnweg+wuphin+wdnhin)*param['dt']
    return(dfs)