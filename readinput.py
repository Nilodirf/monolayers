import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as sp
from scipy import interpolate as ipl
import os
import sys

import florianfit

#this file reads out input parameters from a file and defines additional parameters


def readout():

    class sample:
        def __init__(self, name, gepfit, celfit, cphfit, spin, tc, tdeb, muat, dx, dy, dz, apc, asf, inimag, kappa):
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

    Nickel=sample('Nickel',
                  florianfit.interpol('Ab-initio-Ni/NI_G_ep.txt'),
                  florianfit.interpol('Ab-initio-Ni/NI_c_e.txt'),
                  florianfit.interpol('Ab-initio-Ni/NI_c_p.txt'),
                  0.5, #spin
                  633, #tc
                  360, #tdeb
                  0.616, #muat
                  3.524e-10, 3.524e-10, 3.524e-10, #dxdydz
                  4, #apc
                  0.063,  #asf
                  [0.,0.,0.96925], #inimag at 295K
                  90.7) #kappa
   
    Cobalt=sample('Cobalt',
                  florianfit.interpol('Ab-initio-Co/Co_G_ep.txt'),
                  florianfit.interpol('Ab-initio-Co/Co_c_e.txt'),
                  florianfit.interpol('Ab-initio-Co/Co_c_p.txt'),
                  1.5, #spin
                  1423, #tc
                  342.5, #tdeb
                  1.8, #muat
                  3.54e-10, 3.54e-10, 3.54e-10, #dxdydz
                  4, #apc
                  0.045,  #asf
                  [0.,0.,0.99793], #inimag at 295K
                  90.7)  #kappa

    Iron=sample('Iron',
                florianfit.interpol('Ab-initio-Fe/Fe_G_ep.txt'),
                florianfit.interpol('Ab-initio-Fe/Fe_c_e.txt'),
                florianfit.interpol('Ab-initio-Fe/Fe_c_p.txt'),
                2, #spin
                1041, #tc
                396, #tdeb
                2.2, #muat
                2.856e-10,2.856e-10, 2.856e-10, #dxdydz
                2, #apc
                0.035, #asf
                [0.,0.,0.9863], #inimag at 295K
                0) #kappa

    Gadolinium=sample('Gadolinium',
                      None,
                      None,
                      None,
                      3,
                      293,
                      160,
                      7,
                      3.6e-10,
                      3.6e-10,
                      5.8e-10,
                      6,
                      0.12,
                      [0.,0.,1.],
                      0)

    Dummysd=sample('Dummysd',
                   None,
                   None,
                   None,
                   0.5,
                   700,
                   450,
                   0.7,
                   2.5e-10, 2.5e-10, 2.5e-10,
                   2,
                   None,
                   [0,0,1.],
                   0)
    

    #define different approaches of simulation:
    alexpump=False
    fpflo=False
    model='sd'
    tediff=True
    qes=False

    #sample
    sam=Gadolinium

    #simulation time parameters
    simlen=int(6.4e6)               #length of simulation in units of dt
    dt=1e-16                        #timestep of simulation

    #initial sample conditions
    samplesize=[1,1,1]                  #samplesize
    nj=1                                #film thickness
    h_ext=0                             #external magnetic field
    initemp=100.                        #initial temperature of electron and phonon bath

    #gaussian pulse parameters
    pump_power=1e21                          #power of optical pulse in W/m^3
    pump_sigma=0.0495e-12                    #sigma of gaussian pulse in s
    pump_delay=40e-12                        #position of maximum of gaussian laser pulse in s

    #alternative parameters if no ab initio calculations can be used
    pendep=40e-9                                           #penetration depth of laserpulse
    me0=0.99142842299866007                                #some correction factor for initial magnetization
    cve_const=225                                          #proportionality constant of electron specific heat
    cvp_const=2.33e6#3.78e6                                #constant phonon heat capacity 
    gep=0.25e18                                            #flow constant between electron and phonon heat baths 
    ges=1
    lamda=5e15                                             #constant for ambient heat equilibration of lattice

    #sample constants deducted from class input
    J=float(3*sp.k*sam.tc*(sam.spin**2)/(sam.spin*(sam.spin+1)))                        #exchange splitting constant
    if fpflo:
        R=8*sam.asf*sam.dx*sam.dy*sam.dz/sam.apc*sam.tc**2/sam.tdeb**2/sp.k/sam.muat    #parameter for magnetization dynamics /gep
    else:
        R=8*sam.asf*sam.dx*sam.dy*sam.dz/sam.apc*sam.tc**2/sam.tdeb**2/sp.k/sam.muat
    arbsconst=R/sam.tc**2*J/sp.k                                                    #constant for the computation of arbitrary spin dynamics


    pump=None
    if alexpump:
        pumpfile=open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'3TM_Data/Pump.txt'),'r')
        dpump=np.array([float(i)**2 for i in pumpfile.readlines()])*pump_power
        t=np.arange(len(dpump))*2e-14-2.04e-12+pump_delay
        pump=ipl.interp1d(t,dpump, fill_value=(0,0), bounds_error=False)

    #sample constants for s-d-model
    sdrate=5e12/J*sam.spin
    rhosd=1e-19
    sdissrate=1e13

    #calculate initial kerr-signal:
    #kerr0=0
    #for i in range(samplesize[2]):
    #    kerr0+=math.exp(-(i-1)*sam.dx/pendep)*me0
    #inikerr=kerr0                                               #initial kerr-signal of probe

    param={'gepf': sam.gepfit, 'celf': sam.celfit, 'cphf': sam.cphfit, 'ss': samplesize, 's':sam.spin, 'tc':sam.tc,
           'muat':sam.muat, 'hex':h_ext, 'dt':dt, 'nj':nj, 'R':R, 'inimag':sam.inimag, 'pendep':pendep, 'initemp':initemp, 'me0':me0,
           'dx':sam.dx, 'dy':sam.dy, 'dz':sam.dz, 'simlen':simlen, 'J':J, 'cvec':cve_const, 'cvpc': cvp_const,
           'gep':gep, 'fpflo':fpflo, 'psig':pump_sigma, 'tdeb':sam.tdeb, 'pp':pump_power, 'asc':arbsconst, 'pdel':pump_delay, 
           'asf':sam.asf, 'mod':model, 'apc':sam.apc, 'name':sam.name, 'model':model, 'pump':pump, 'ap':alexpump, 'kappa':sam.kappa, 'tediff':tediff, 'ges':ges, 'qes':qes,
           'lambda':lamda, 'sdrate':sdrate, 'rhosd':rhosd, 'sdissrate':sdissrate}

    if nj<samplesize[2]:
        print('The magnetic sampledepth (nz) must not be larger then the total (thermally excited) sampledepth (nj)')
    return(param)

param=readout()
