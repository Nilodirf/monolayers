import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as sp

from readinput import param


def tempdyn(ts, te, tp, dqes):

    if param['fpflo']:
        cv_el=param['celf'](te)
        cv_ph=param['cphf'](tp)
    else:
        #compute electronic specific heat:#
        cv_el=te*param['cvec']
        cv_ph=param['cvpc']
        if param['name']=='Gadolinium':
            cv_ph=(1.51e6*(120/tp)**2)*np.exp(120/tp)/((np.exp(120/tp)-1)**2) 
        #compute phononic specific heat depending on temperature#
        #cv_ph=np.zeros(param['nj'])
        #tred=np.true_divide(param['tdeb'],tp)
        #trede=np.true_divide(param['tein'],tp)

        ##define boolean arrays to evaluate only necessary entries:
        #masktred=tred>3.5                                       
        #nomasktred=np.logical_not(masktred) 
    
        ##Taylor expansion:
        #cv_ph[masktred]=np.sum(tp[masktred,np.newaxis]**np.arange(1,param['dp'].size+1)*param['dp'],axis=-1) 
    
        ##Debye integral:
        #cv_ph[masktred]=debint(tred[masktred])   

        ##evaluate other entires not fulfilling condidion masktred                                                                  
        #cv_ph[np.logical_not(masktred)]=2.55e6*np.power(trede[nomasktred],2)*np.true_divide(np.exp(trede[nomasktred]),np.power(np.exp(trede[nomasktred])-1,2))

        ##according to Florians simulation:
        #cv_ph=param['cvpc']


    ############add magnon specific heat if T<Tc#############
 
    #cv_mag=np.zeros(param['nj'])
    #masktc = tp<param['tc']
    #cv_mag[masktc] = np.sum(param['mp']*(tp[masktc,np.newaxis]**np.arange(param['mp'].size)),axis=-1)
    #cv_ph+=cv_mag

    ###########compute pump temperature for each timestep and depth:##########
    if param['ap']:
        dpump=param['pump'](ts*param['dt'])
    else:
        dpump=param['pp']*np.exp(-np.arange(param['nj'])*param['dz']/param['pendep'])*np.exp(-((ts*param['dt']-param['pdel'])**2)/2/param['psig']**2)

    ######compute heat diffusion#####
        #if ts>9e4 and ts%1000==0:
        #    print(tediff*param['kappa']*param['dz']**2)

    ############calculate temperature derivatives##########

    if param['fpflo']:
        teq=param['gepf'](te)*(te-tp)
    else:
        teq=param['gep']*(te-tp)

    dtel=dpump-teq
    dtph=teq
    if param['qes']:
        dtel+=+param['ges']*dqes
        dtph+=(1-param['ges'])*dqes
    if param['tediff']:
        if param['name']=='Gadolinium':
            #dtel=dtel+param['kappa']*te/tp*tediff/param['dz']**2
            dtph+=param['lambda']*(100-tp)
        else:
            tenext=np.roll(te,1)
            telast=np.roll(te,-1)
            tenext[0]=0
            telast[-1]=0
            tediff=tenext+telast-2*te
            tediff[0]+=te[0]
            tediff[-1]+=te[-1]
            dtel=dtel+param['kappa']*tediff/param['dz']**2
    dte=param['dt']*np.true_divide(dtel,cv_el)
    dtp=param['dt']*np.true_divide(dtph,cv_ph)

    return(np.array(dte), np.array(dtp), dpump)
    

#computation of the debye integral:

def debint(b):
    grains=2**8
    delta=np.true_divide(b,grains)
    x=delta
    sum=np.zeros(param['nj'],dtype=np.float64)
    for gr in range(1,grains):
        sum=sum+np.multiply(integrand(x),delta)
        x+=delta
    cp=0.5*3*2.55e6*((b)**(-3))+sum
    return(cp)

def integrand(y):
    num=np.power(y,4)*np.exp(y)
    rdenom=np.exp(y)-1
    denom=rdenom**2
    quot=np.true_divide(num,denom)  
    return(quot)