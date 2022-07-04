import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as sp
import os
import sys
from datetime import datetime

#import other files
import magdyn
import tempdyn
import magdyn_s
import sd_mag
import pf_mag
date=datetime.now().replace(microsecond=0)

#this file defines the initial parameters in a usable form and calls the function(s) to compute the dynamics

def output(pl, pp):
    #create file and document input data
    file=open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), '3TM_results/'+ str(pl['name']) +'/finite depth/' + str(pp) + '.dat'), 'w+')

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
        file.write('# gep=' +str(pl['gepf']) + '[W/m^3/K]' + '\n')
        file.write('# factor cv_el=' + str(pl['celf']) + '[J/m^3/K^2]' + '\n')
        if pl['name']=='Gadolinium':
            file.write('# cv_ph computed with Einstein model, T_Ein=0.75*' + str(pl['tdeb']) + '[K]' + '\n')
        else:
            file.write('# cv_ph=' + str(pl['cphf']) + '[J/m^3/K]' + '\n')
        file.write('# R/g_ep = ' + str(pl['R']) + '[1/s]' + '\n')
    file.write('# Tc=' + str(pl['tc']) + '[K]' + '\n')
    file.write('# T_Deb=' + str(pl['tdeb']) + '[K]' + '\n')
    file.write('### pulse paramters ###' + '\n')
    file.write('# peak power:' + str(pl['pp']) + '[W/m^3]' + '\n')
    file.write('# sigma=' + str(pl['psig']) + '[s]' + '\n')
    file.write('# delay=' + str(pl['pdel']) + '[s]' + '\n')
    file.write('#############################' + '\n')
    if pl['model']=='sd':
        file.write('# time [ps]'+'\t'+'magnetization'+'\t'+'T_e [K]'+'\t'+'T_p [K]'+'\t'+'E_tot [J/m^3]'+'\t'
                   +'E_el [J/m^3]'+'\t'+'E_ph [J/m^3]'+'\t'+'E_spin [J/m^3]'+ '\t' + 'mu_s [J]' + '\n')

    elif pl['model']=='pf':
        file.write('# time [ps]' + '\t' + 'magnetization' + '\t' + 'T_e [K]' + '\t' + 'T_p [K]' + '\t' + 'E_tot [J/m^3]' + '\t'
                    + 'E_el [J/m^3]' + '\t' + 'E_ph [J/m^3]' + '\t' + 'E_spin [J/m^3]' + '\t' + 'it. magnetization' + '\n')

    else:
        file.write('# time [ps]'+'\t'+'magnetization'+'\t'+'T_e [K]'+'\t'+'T_p [K]'+'\t'+'E_tot [J/m^3]'+'\t'
                +'E_el [J/m^3]'+'\t'+'E_ph [J/m^3]'+ '\t' +'E_spin [J/m^3]'+ '\n')

    ###initiate stuff for any model###
    samsize=[pl['ss'][i] for i in range(3)]
    lattice=np.arange(samsize[2])
    exp0=np.array([np.exp(-i*pl['dz']/pl['pendep']) for i in lattice])
    dqes = np.array([0.])
    
    ###initialize absorbed energy of each subsystem###
    etot = 0
    eel = 0
    eph = 0
    espin = 0

    tempe = np.array([pl['initemp'] for i in range(pl['nj'])])
    tempp = np.array([pl['initemp'] for i in range(pl['nj'])])


    if pl['model']=='m3tm':

        magx=np.array([[[pl['inimag'][0] for i in range(samsize[2])] for j in range(samsize[1])] for k in range(samsize[0])])
        magy=np.array([[[pl['inimag'][1] for i in range(samsize[2])] for j in range(samsize[1])] for k in range(samsize[0])])
        magz=np.array([pl['inimag'][2] for i in range(samsize[2])])

        kerr0 = np.sum(magz * exp0)

        for t in range(pl['simlen']):

            if t==0 or t%100==0:
                kerr = np.sum(magz * exp0) / kerr0
                file.write(str(t * pl['dt'] * 1e12) + '\t' + str(kerr) + '\t' + str(tempe[0]) + '\t' + str(tempp[0]) + '\t'
                + str(etot) + '\t' + str(eel) + '\t' + str(eph) + '\t' + str(espin) + '\n')

            dmagz = magdyn.magdyn(t, samsize, magz, tempe, tempp)
            magzerr = magz + dmagz / 2
            magz2 = magz + dmagz
            dmagz2 = magdyn.magdyn(t, samsize, magz2, tempe, tempp)
            magz = magzerr + dmagz2 / 2

            if pl['qes'] and t > pl['pdel'] / pl['dt'] - 1e4:
                dqes = (dmagz + dmagz2) / 2 / pl['dt'] * pl['J'] * magz / pl['dx'] / pl['dy'] / pl['dz'] * pl['apc']

            dtemps = tempdyn.tempdyn(t, tempe, tempp, dqes, pp)
            tempe += dtemps[0]
            tempp += dtemps[1]

            eel += pl['celf'](tempe[0]) * dtemps[0][0]
            deel = pl['celf'](tempe[0]) * dtemps[0][0]
            eph += pl['cphf'](tempp[0]) * dtemps[1][0]
            deph = pl['cphf'](tempp[0]) * dtemps[1][0]
            espin += dqes[0] * pl['dt']
            despin = dqes[0] * pl['dt']
            etot += deel + deph - despin

        file.close()

    elif pl['model']=='arbspin':

        fs0 = np.zeros(int(2 * pl['s']) + 1)
        fs0[0] = 1
        # for an initial magnetization of 1, all other states ms are unoccupied:
        fs = np.array([fs0 for i in range(samsize[2])])
        print(fs)
        ms = (np.arange(2 * pl['s'] + 1) + np.array([-pl['s'] for i in range(int(2 * pl['s']) + 1)]))
        print(ms)
        sup = -np.power(ms, 2) - ms + pl['s'] ** 2 + pl['s']
        print(sup)
        sdn = -np.power(ms, 2) + ms + pl['s'] ** 2 + pl['s']
        print(sdn)
        magz = -np.sum(ms * fs, axis=-1) / pl['s']

        kerr0 = np.sum(magz * exp0)

        for t in range(pl['simlen']):

            if t==0 or t%100==0:
                kerr = np.sum(magz * exp0) / kerr0
                file.write(str(t * pl['dt'] * 1e12) + '\t' + str(kerr) + '\t' + str(tempe[0]) + '\t' + str(tempp[0]) + '\t'
                + str(etot) + '\t' + str(eel) + '\t' + str(eph) + '\t' + str(espin) + '\n')

            dfs = magdyn_s.magdyn_s(t, samsize, magz, fs, sup, sdn, tempe, tempp)
            dmagz = -np.sum(ms * dfs, axis=-1) / pl['s']
            fs2 = fs + dfs
            magz2 = magz + dmagz
            dfs2 = magdyn_s.magdyn_s(t, samsize, magz2, fs2, sup, sdn, tempe, tempp)
            dmagz2 = -np.sum(ms * dfs2, axis=-1) / pl['s']
            df = (dfs + dfs2) / 2
            dm = (dmagz + dmagz2) / 2
            # if t*pl['dt']-pl['pdel']<=100e-15:
            #    dm=2*dm
            #    df=2*df
            newfs = fs + df
            newmagz = magz + dm
            fs = newfs
            magz = newmagz

            if pl['qes'] and t > pl['pdel'] / pl['dt'] - 1e4:
                dqes = (dmagz + dmagz2) / 2 / pl['dt'] * pl['J'] * magz / pl['dx'] / pl['dy'] / pl['dz'] * pl['apc']

            dtemps = tempdyn.tempdyn(t, tempe, tempp, dqes)
            tempe += dtemps[0]
            tempp += dtemps[1]

            eel += pl['celf'](tempe[0]) * dtemps[0][0]
            deel = pl['celf'](tempe[0]) * dtemps[0][0]
            eph += pl['cphf'](tempp[0]) * dtemps[1][0]
            deph = pl['cphf'](tempp[0]) * dtemps[1][0]
            espin += dqes[0] * pl['dt']
            despin = dqes[0] * pl['dt']
            etot += deel + deph - despin

        file.close()

    elif pl['model']=='sd':

        fs0 = np.zeros(int(2 * pl['s']) + 1)
        fs0[0] = 1
        fs = np.array([fs0 for i in range(samsize[2])])
        ms = (np.arange(2 * pl['s'] + 1) + np.array([-pl['s'] for i in range(int(2 * pl['s']) + 1)]))
        sup = -np.power(ms, 2) - ms + pl['s'] ** 2 + pl['s']
        sdn = -np.power(ms, 2) + ms + pl['s'] ** 2 + pl['s']
        magz = -np.sum(ms * fs, axis=-1) / pl['s']
        mus = np.array([0 for i in range(samsize[2])])

        kerr0 = np.sum(magz * exp0)

        for t in range(pl['simlen']):

            if t==0 or t%100==0:
                kerr = np.sum(magz * exp0) / kerr0
                file.write(str(t * pl['dt'] * 1e12) + '\t' + str(kerr) + '\t' + str(tempe[0]) + '\t' + str(tempp[0]) + '\t'
                + str(etot) + '\t' + str(eel) + '\t' + str(eph) + '\t' + str(espin) + str(mus[0]) + '\n')

            dfs = sd_mag.locmag(samsize, magz, tempe, mus, fs, sup, sdn, pf, pp)
            dmagz = -np.sum(ms * dfs, axis=-1) / pl['s']
            magzz2 = magz + dmagz
            dmus = sd_mag.itmag(samsize, tempe, dmagz, mus, pf, pp)
            mus2 = mus + dmus
            dfs2 = sd_mag.locmag(samsize, magzz2, tempe, mus2, fs, sup, sdn, pf, pp)
            newfs = fs + (dfs + dfs2) / 2
            dmagz2 = -np.sum(ms * dfs2, axis=-1) / pl['s']
            newmagz = magz + (dmagz + dmagz2) / 2
            dmus2 = sd_mag.itmag(samsize, tempe, dmagz2, mus2, pf, pp)
            newmus = mus + (dmus + dmus2) / 2
            magz = newmagz
            fs = newfs
            mus = newmus

            if pl['qes'] and t > pl['pdel'] / pl['dt'] - 1e4:
                dqes = (dmagz + dmagz2) / 2 / pl['dt'] * pl['J'] * magz / pl['dx'] / pl['dy'] / pl['dz'] * pl['apc']

            dtemps = tempdyn.tempdyn(t, tempe, tempp, dqes)
            tempe += dtemps[0]
            tempp += dtemps[1]

            eel += pl['celf'](tempe[0]) * dtemps[0][0]
            deel = pl['celf'](tempe[0]) * dtemps[0][0]
            eph += pl['cphf'](tempp[0]) * dtemps[1][0]
            deph = pl['cphf'](tempp[0]) * dtemps[1][0]
            espin += dqes[0] * pl['dt']
            despin = dqes[0] * pl['dt']
            etot += deel + deph - despin

        file.clos()

    elif pl['model'] == 'pf':

        fs0 = np.zeros(int(2 * pl['s']) + 1)
        fs0[0] = 1
        fs = np.array([fs0 for i in range(samsize[2])])
        ms = (np.arange(2 * pl['s'] + 1) + np.array([-pl['s'] for i in range(int(2 * pl['s']) + 1)]))
        sup = -np.power(ms, 2) - ms + pl['s'] ** 2 + pl['s']
        sdn = -np.power(ms, 2) + ms + pl['s'] ** 2 + pl['s']
        magz = -np.sum(ms * fs, axis=-1) / pl['s']
        magp = np.array([0 for i in range(samsize[2])])

        kerr0 = np.sum(magz * exp0)

        for t in range(pl['simlen']):

            if t == 0 or t % 100 == 0:
                kerr = np.sum(magz * exp0) / kerr0
                kerrp = np.sum(magp*exp0) / kerr0
                file.write(
                    str(t * pl['dt'] * 1e12) + '\t' + str(kerr) + '\t' + str(tempe[0]) + '\t' + str(tempp[0]) + '\t'
                    + str(etot) + '\t' + str(eel) + '\t' + str(eph) + '\t' + str(espin) + '\t' + str(kerrp) + '\n')

            dfs = pf_mag.locmag(samsize, magz, tempe, magp, fs, sup, sdn, pf, pp)
            dmagz = -np.sum(ms * dfs, axis=-1) / pl['s']
            magzz2 = magz + dmagz
            dmagp = pf_mag.itmag(samsize, tempe, dmagz, magp, pf, pp)
            magp2 = magp + dmagp
            dfs2 = pf_mag.locmag(samsize, magzz2, tempe, magp2, fs, sup, sdn, pf, pp)
            newfs = fs + (dfs + dfs2) / 2
            dmagz2 = -np.sum(ms * dfs2, axis=-1) / pl['s']
            newmagz = magz + (dmagz + dmagz2) / 2
            dmagp2 = pf_mag.itmag(samsize, tempe, dmagz2, magp2, pf, pp)
            newmagp = magp + (dmagp + dmagp2) / 2
            magz = newmagz
            fs = newfs
            magp = newmagp

            if pl['qes'] and t > pl['pdel'] / pl['dt'] - 1e4:
                dqes = (dmagz + dmagz2) / 2 / pl['dt'] * pl['J'] * magz / pl['dx'] / pl['dy'] / pl['dz'] * pl['apc']

            dtemps = tempdyn.tempdyn(t, tempe, tempp, dqes)
            tempe += dtemps[0]
            tempp += dtemps[1]

            eel += pl['celf'](tempe[0]) * dtemps[0][0]
            deel = pl['celf'](tempe[0]) * dtemps[0][0]
            eph += pl['cphf'](tempp[0]) * dtemps[1][0]
            deph = pl['cphf'](tempp[0]) * dtemps[1][0]
            espin += dqes[0] * pl['dt']
            despin = dqes[0] * pl['dt']
            etot += deel + deph - despin

        file.close()
        
    return