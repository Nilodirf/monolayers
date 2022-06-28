import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as sp
from scipy import interpolate as ipl
import os
import sys

import florianfit


def readout():

    class sample:
        def __init__(self, name, gepfit, celfit, cphfit, spin, tc, tdeb, muat, dx, dy, dz, apc, asf, inimag, kappa, locmom, locspin, model):
            self.name=name
            self.gepfit=gepfit
            self.celfit=celfit
            self.cphfit=cphfit
            self.spin=spin
            self.tc=tc
            self.tdeb=tdeb
            self.muat=muat
            self.dx=dx
            self.dy=dy
            self.dz=dz
            self.apc=apc
            self.asf=asf
            self.inimag=inimag
            self.kappa=kappa
            self.locmom=locmom
            self.locspin=locspin
            self.model=model
            self.J=float(3*sp.k*self.tc*((self.spin-self.locspin)**2)/((self.spin-self.locspin)*((self.spin-self.locspin)+1)))
            self.Jloc = float(3 * sp.k * self.tc * (self.locspin ** 2) / (self.locspin * (self.locspin + 1)))
            self.R =8*sam.asf*sam.dx*sam.dy*sam.dz/sam.apc*sam.tc**2/sam.tdeb**2/sp.k/(sam.muat-sam.locmom)
            self.arbsc=self.R/sam.tc**2*self.J/sp.k


    Nickel=sample('Nickel',
                  florianfit.interpol('Ab-initio-Ni/Ni_G_ep.txt', None, None),
                  florianfit.interpol('Ab-initio-Ni/Ni_c_e.txt', None, None),
                  florianfit.interpol('Ab-initio-Ni/Ni_c_p.txt', None, None),
                  0.5, #spin
                  633, #tc
                  360, #tdeb
                  0.616, #muat
                  3.524e-10, 3.524e-10, 3.524e-10, #dxdydz
                  4, #apc
                  0.063,  #asf
                  [0.,0.,0.96925], #inimag at 295K
                  90.7, #kappa
                  0.,
                  0.,
                  'm3tm')
   
    Cobalt=sample('Cobalt',
                  florianfit.interpol('Ab-initio-Co/Co_G_ep.txt', None, None),
                  florianfit.interpol('Ab-initio-Co/Co_c_e.txt', None, None),
                  florianfit.interpol('Ab-initio-Co/Co_c_p.txt', None, None),
                  1.5, #spin
                  1423, #tc
                  342.5, #tdeb
                  1.8, #muat
                  3.54e-10, 3.54e-10, 3.54e-10, #dxdydz
                  4, #apc
                  0.045,  #asf
                  [0.,0.,0.99793], #inimag at 295K
                  90.7, #kappa
                  0.,
                  0.,
                  'arbspin')

    Iron=sample('Iron',
                florianfit.interpol('Ab-initio-Fe/Fe_G_ep.txt', None, None),
                florianfit.interpol('Ab-initio-Fe/Fe_c_e.txt', None, None),
                florianfit.interpol('Ab-initio-Fe/Fe_c_p.txt', None, None),
                2, #spin
                1041, #tc
                396, #tdeb
                2.2, #muat
                2.856e-10,2.856e-10, 2.856e-10, #dxdydz
                2, #apc
                0.035, #asf
                [0.,0.,0.9863], #inimag at 295K
                0, #kappa
                0.,
                0.,
                'arbspin')

    Gadolinium=sample('Gadolinium',
                      florianfit.interpol(2.5e17, 'const', None),
                      florianfit.interpol(225., 'lin', None),
                      florianfit.interpol(1.51e6, 'einstein', 120),
                      3.5,
                      293,
                      160,
                      7.5,
                      3.6e-10,
                      3.6e-10,
                      5.8e-10,
                      7.5,
                      0.12,
                      [0.,0.,1.],
                      0,
                      7.,
                      3.,
                      'sd')
    

    #define different approaches of simulation:
    alexpump=False
    tediff=False
    fpflo=False
    qes=False

    #sample
    sam=Gadolinium

    #simulation time parameters
    simlen=int(1e6)               #length of simulation in units of dt
    dt=1e-16                        #timestep of simulation

    #initial sample conditions
    samplesize=[1,1,1]                  #samplesize
    nj=1                                #film thickness
    h_ext=0                             #external magnetic field
    initemp=100.                        #initial temperature of electron and phonon bath [K]

    #gaussian pulse parameters
    pump_power=1e21                          #power of optical pulse in W/m^3
    pump_sigma=0.0495e-12                    #sigma of gaussian pulse in s
    pump_delay=10e-12                        #position of maximum of gaussian laser pulse in s
    pendep=40e-9                             #penetration depth of laserpulse

    ges=1
    lamda=5e15                                             #constant for ambient heat equilibration of lattice

    pump=None
    if alexpump:
        pumpfile=open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'3TM_Data/Pump.txt'),'r')
        dpump=np.array([float(i)**2 for i in pumpfile.readlines()])*pump_power
        t=np.arange(len(dpump))*2e-14-2.04e-12+pump_delay
        pump=ipl.interp1d(t,dpump, fill_value=(0,0), bounds_error=False)

    #sample constants for s-d-model
    sdrate=np.array([0.1])*1e-12
    rhosd=1e-19
    sdissrate=np.array([0.1])*1e-12

    #calculate initial kerr-signal:
    #kerr0=0
    #for i in range(samplesize[2]):
    #    kerr0+=math.exp(-(i-1)*sam.dx/pendep)*me0
    #inikerr=kerr0                                               #initial kerr-signal of probe

    param={'gepf': sam.gepfit, 'celf': sam.celfit, 'cphf': sam.cphfit, 'ss': samplesize, 's':sam.spin, 'tc':sam.tc,
           'muat':sam.muat, 'hex':h_ext, 'dt':dt, 'nj':nj, 'R':sam.R, 'inimag':sam.inimag, 'pendep':pendep, 'initemp':initemp,
           'dx':sam.dx, 'dy':sam.dy, 'dz':sam.dz, 'simlen':simlen, 'J':sam.J,'fpflo':fpflo, 'psig':pump_sigma, 'tdeb':sam.tdeb, 'pp':pump_power,
           'asc':sam.arbsc, 'pdel':pump_delay,
           'asf':sam.asf, 'apc':sam.apc, 'name':sam.name, 'model':sam.model, 'pump':pump, 'ap':alexpump, 'kappa':sam.kappa, 'tediff':tediff, 'ges':ges, 'qes':qes,
           'lambda':lamda, 'sdrate':sdrate, 'rhosd':rhosd, 'sdissrate':sdissrate, 'locspin':sam.locspin, 'locmom':sam.locmom}

    if nj<samplesize[2]:
        print('The magnetic sampledepth (nz) must not be larger then the total (thermally excited) sampledepth (nj)')
    return(param)

param=readout()
