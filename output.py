import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as sp
from datetime import datetime

#import other files
import magdyn
import tempdyn
import magdyn_s
import sd_mag
date=datetime.now().replace(microsecond=0)

#this file defines the initial parameters in a usable form and calls the function(s) to compute the dynamics

def output(pl):

    #create file and document input data
    file=open('C:/Users/tgrie/Desktop/3TM_results/'+str(pl['name'])+'/'+str(input())+'.dat','w+')

    file.write('# time of execution: ' + str(date) + '\n')
    file.write('#model:' + str(pl['model']) + '\n')
    file.write('### simulation parameters ###' + '\n')
    file.write('# timestep dt=' + str(pl['dt']) + '[s]' + '\n')
    file.write('# initial temperature:' + str(pl['initemp']) + '[K]' + '\n')
    file.write('# initial magnetization [x,y,z]:' + str(pl['inimag']) + '\n')
    file.write('### sample paramters ###' + '\n')
    file.write('# Sample: ' + str(pl['name']) + '\n')
    file.write('# S=' + str(pl['s']) + '\n')
    file.write('# mu_at=' + str(pl['muat']) + '[mu_B]' + '\n')
    file.write('# nx=' + str(pl['ss'][0]) + '\t' + 'ny=' + str(pl['ss'][1]) + '\t' + 'nz=' + str(pl['ss'][2]) + '\n')
    file.write('# dx=' + str(pl['dx']) + '\t' + 'dy=' + str(pl['dy']) + '\t' + 'dz=' + str(pl['dz']) + '\t' + '[m]' '\n')
    file.write('# atoms per unit cell:' + str(pl['apc']) + '\n')
    file.write('# asf=' + str(pl['asf']) + '\n')
    if pl['tediff']:
        file.write('# lambda=' + str(pl['lambda']) + '[W/m^3/K]' + '\n')
    if pl['fpflo']:
        file.write('# gep, cv_el, cv_ph: DFT-interpolation' + '\n')
    else:
        file.write('# gep=' +str(pl['gep']) + '[W/m^3/K]' + '\n')
        file.write('# factor cv_el=' + str(pl['cvec']) + '[J/m^3/K^2]' + '\n')
        if pl['name']=='Gadolinium':
            file.write('# cv_ph computed with Einstein model, T_Ein=0.75*' + str(pl['tdeb']) + '[K]' + '\n')
        else:
            file.write('# cv_ph=' + str(pl['cvpc']) + '[J/m^3/K]' + '\n')
        file.write('# R/g_ep = ' + str(pl['R']) + '[1/s]' + '\n')
    file.write('# Tc=' + str(pl['tc']) + '[K]' + '\n')
    file.write('# T_Deb=' + str(pl['tdeb']) + '[K]' + '\n')
    file.write('### pulse paramters ###' + '\n')
    file.write('# peak power:' + str(pl['pp']) + '[W/m^3]' + '\n')
    file.write('# sigma=' + str(pl['psig']) + '[s]' + '\n')
    file.write('# delay=' + str(pl['pdel']) + '[s]' + '\n')
    file.write('#############################' + '\n')
    file.write('# time [ps]'+'\t'+'magnetization'+'\t'+'T_e [K]'+'\t'+'T_p [K]'+'\t'+'E_tot [J/m^3]'+'\t'
               +'E_el [J/m^3]'+'\t'+'E_ph [J/m^3]'+'\t'+'E_spin [J/m^3]'+ '\n')


    samsize=[pl['ss'][i] for i in range(3)]
    lattice=np.arange(samsize[2])
    exp0=np.array([np.exp(-i*pl['dz']/pl['pendep']) for i in lattice])
    
    #####magnetization initialization######
    if pl['model']=='m3tm':
        magx=np.array([[[pl['inimag'][0] for i in range(samsize[2])] for j in range(samsize[1])] for k in range(samsize[0])])
        magy=np.array([[[pl['inimag'][1] for i in range(samsize[2])] for j in range(samsize[1])] for k in range(samsize[0])]) 
        magz=np.array([pl['inimag'][2] for i in range(samsize[2])])
    #for arbitrary spin:
    elif pl['model']=='arbspin':
        fs0=np.zeros(int(2*pl['s'])+1)
        fs0[0]=1
        #for an initial magnetization of 1, all other states ms are unoccupied:
        fs=np.array([fs0 for i in range(samsize[2])])
        print(fs)
        ms=(np.arange(2*pl['s']+1)+np.array([-pl['s'] for i in range(int(2*pl['s'])+1)]))
        print(ms)
        sup=-np.power(ms,2)-ms+pl['s']**2+pl['s']
        print(sup)
        sdn=-np.power(ms,2)+ms+pl['s']**2+pl['s']
        print(sdn)
        magz=-np.sum(ms*fs, axis=-1)/pl['s']
    elif pl['model']=='sd':
        #fs0=np.zeros(int(2*pl['s'])+1)
        #fs0[0]=1
        #fs=np.array([[[fs0 for i in range(samsize[2])] for j in range(samsize[1])] for k in range(samsize[0])])
        #ms=(np.arange(2*pl['s']+1)+np.array([-pl['s'] for i in range(int(2*pl['s'])+1)]))
        #sup=-np.power(ms,2)-ms+pl['s']**2+pl['s']
        #sdn=-np.power(ms,2)+ms+pl['s']**2+pl['s']
        #magz=-np.sum(ms*fs, axis=-1)/pl['s']
        magz=np.array([pl['inimag'][2] for i in range(samsize[2])])
        mus=np.array([0 for i in range(samsize[2])])
    dqes=np.array([0.])
    kerr0=np.sum(magz*exp0)
    ###initialize absorbed energy of each subsystem###
    etot=0
    eel=0
    eph=0
    espin=0

    dm=np.array([0])


    ####temperature initialisation:#####
    tempe=np.array([pl['initemp'] for i in range(pl['nj'])])
    tempp=np.array([pl['initemp'] for i in range(pl['nj'])])

    
    for t in range(pl['simlen']):

        kerr=np.sum(magz*exp0)/kerr0

        if t==0 or t%100==0:
            file.write(str(t*pl['dt']*1e12) + '\t' + str(kerr) + '\t' + str(tempe[0]) + '\t' + str(tempp[0]) + '\t' 
                     + str(etot) + '\t' + str(eel) + '\t' + str(eph) + '\t' + str(espin) + '\n')

        #if t==0 or t%100==0:
        #    file.write(str(t*pl['dt']*1e12) + '\t' + str(magz[0,0,0]/pl['inimag'][2]) + '\t' + str(tempe[0]) + '\t' + str(tempp[0]) + '\n')

        ######this block computes the magnetization dynamics in the Heun method######
        if pl['model']=='m3tm':
            dmagz=magdyn.magdyn(t, samsize, magz, tempe, tempp)
            magzerr=magz+dmagz/2
            magz2=magz+dmagz
            dmagz2=magdyn.magdyn(t, samsize, magz2, tempe, tempp)
            magz=magzerr+dmagz2/2
            #if t%100==0:
                #dm=0
            #else:
                #dm+=dmagz+dmagz2
        ####this block computes the magnetization dynamics for arbitrary spin (with Euler method)####
        elif pl['model']=='arbspin':
            dfs=magdyn_s.magdyn_s(t, samsize, magz, fs, sup, sdn, tempe, tempp)
            dmagz=-np.sum(ms*dfs, axis=-1)/pl['s']
            fs2=fs+dfs
            magz2=magz+dmagz
            dfs2=magdyn_s.magdyn_s(t, samsize, magz2, fs2, sup, sdn, tempe, tempp)
            
            dmagz2=-np.sum(ms*dfs2, axis=-1)/pl['s']
            df=(dfs+dfs2)/2
            dm=(dmagz+dmagz2)/2
            #if t*pl['dt']-pl['pdel']<=100e-15:
            #    dm=2*dm
            #    df=2*df
            newfs=fs+df
            newmagz=magz+dm
            fs=newfs
            magz=newmagz
        elif pl['model']=='sd':
            #dfs=sd_mag.sd_mag(t, samsize, magz, tempe, mus)
            #fserr=fs+dfs/2
            #dmagz=-np.sum(ms*dfs, axis=-1)/pl['s']
            #magzerr=magz+dmagz/2
            #fs2=fs+dfs
            #magz2=magz+dmagz
            #dfs2=sd_mag.sd_mag(t, samsize, magz2, tempe, sup, sdn)
            #newfs=fserr+dfs2/2
            #dmagz2=-np.sum(ms*dfs2, axis=-1)/pl['s']
            #newmagz=magzerr+dmagz2/2
            dmagz=sd_mag.sd_mag(t, samsize, magz, tempe, mus)[0]
            magzerr=magz+dmagz/2
            magz2=magz+dmagz
            dmus=sd_mag.sd_mag(t, samsize, magz, tempe, mus)[1]
            muserr=mus+dmus/2
            mus2=mus+dmus
            dmagz2=sd_mag.sd_mag(t, samsize, magz2, tempe, mus2)[0]
            magz=magzerr+dmagz2/2
            dmus2=sd_mag.sd_mag(t, samsize, magz2, tempe, mus2)[1]
            mus=muserr+dmus2/2
        
        if pl['qes'] and t>pl['pdel']/pl['dt']-1e4:
            dqes=(dmagz+dmagz2)/2/pl['dt']*pl['J']*magz/pl['dx']/pl['dy']/pl['dz']*pl['apc']
            #x=magz*pl['J']/sp.k/pl['tc']
            #heff=sp.k*tempe[:samsize[2]]/pl['muat']*(1-magdyn.brillouin(x,pl['s'])/magz)/magdyn.dbrillouin(x,pl['s'])*magz
            #dqes=(dmagz+dmagz2)/2/pl['dt']*heff
        dtemps=tempdyn.tempdyn(t, tempe, tempp, dqes)
        tempe+=dtemps[0]
        tempp+=dtemps[1]

        #####this integrates the absorbed power in each subsystem####f
        if pl['fpflo']:
            eel+=pl['celf'](tempe[0])*dtemps[0][0]
            deel=pl['celf'](tempe[0])*dtemps[0][0]
            eph+=pl['cphf'](tempp[0])*dtemps[1][0]
            deph=pl['cphf'](tempp[0])*dtemps[1][0]
        else:
            eel+=pl['cvec']*tempe[0]*dtemps[0][0]
            deel=pl['cvec']*tempe[0]*dtemps[0][0]
            eph+=pl['cvpc']*dtemps[1][0]
            deph=pl['cvpc']*dtemps[1][0]
        espin+=dqes[0]*pl['dt']
        despin=dqes[0]*pl['dt']
        etot+=deel+deph-despin
    file.close()
    return