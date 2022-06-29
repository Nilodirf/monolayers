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
    itmom=param['muat']-param['locmom']
    p_coup= itmom*param['J']*mz/(param['s']-param['locspin'])
    f_coup= param['locmom']*param['Jloc']*mp/param['locspin']
    const_ijk=1/pf*(f_coup-p_coup)/sp.k/param['tc']/param['locmom']/np.sinh((f_coup-p_coup)/param['locmom']/(2*sp.k*param['tc']))
    fsup=sup*fs
    fsdn=sdn*fs
    wupwegprep=const_ijk*np.exp(-(f_coup-p_coup)/param['locmom']/(2*sp.k*te[:ss[2]]))
    wdnwegprep=const_ijk*np.exp((f_coup-p_coup)/param['locmom']/(2*sp.k*te[:ss[2]]))
    wupweg=wupwegprep[...,np.newaxis]*fsup
    wdnweg=wdnwegprep[...,np.newaxis]*fsdn
    wuphin=np.roll(wupweg,1)
    wdnhin=np.roll(wdnweg,-1)
    dfs=(-wupweg-wdnweg+wuphin+wdnhin)*param['dt']
    return(dfs)