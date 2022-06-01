import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import constants as sp

from readinput import param


#this file computes the magnetization daynamics for each time step

def magdyn(ts, ss, mz, te, tp):
        #compute M3TM magnetization dynamics for z component of magnetization
        hmfaz=param['J']*mz+param['hex']
        eta=hmfaz/sp.k/te[:ss[2]]
        mmag=brillouin(eta,param['s'])
        #mmag=np.tanh(param['tc']/te[:ss[2]]*mz)
        if param['fpflo']:
            qs=3*param['tc']*mmag/2/(param['s']+1)/te[:ss[2]]
            gamma=tp[:ss[2]]/param['tc']*param['R']*param['gepf'](te[:ss[2]])#*1/(2*qs/np.sinh(2*qs))
        else:
            #qs=3*param['tc']*mmag/2/(param['s']+1)/te[:ss[2]]
            gamma=tp[:ss[2]]/param['tc']*param['R']*param['gep']#*1/(2*qs/np.sinh(2*qs))*1.3
        rate=(1-np.true_divide(mz,mmag))*gamma
        dmz=rate*mz*param['dt']
        return(dmz)




def brillouin(x,spin):
    #Compute equilibrium magnetization via Brillouin function
    c1=(2*spin+1)/(2*spin)
    c2=1/(2*spin)
    fb=c1/np.tanh(c1*x)-c2/np.tanh(c2*x)
    return(fb)

def dbrillouin(x,spin):
    c1=(2*spin+1)/(2*spin)
    c2=1/(2*spin)
    dfb=c2**2/(np.sinh(c2*x))**2-c1**2/(np.sinh(c1*x))**2
    return(dfb)