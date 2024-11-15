#!/usr/bin/env python3

import numpy as np
from scipy.fft import fft
import os
import sys, getopt
from progress.bar import Bar

inputfile = ''
outputfile = ''
try:
    
    opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["ifile=","ofile="])
    if len(opts) != 2:
        print ('Spectra_1.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
except getopt.GetoptError:
    print ('Spectra_1.py -i <inputfile> -o <outputfile>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print ('Spectra_1.py -i <inputfile> -o <outputfile>')
        sys.exit()
    elif opt in ("-i", "--ifile"):
        inputfile = arg
    elif opt in ("-o", "--ofile"):
        outputfile = arg
   



Data=np.load(inputfile)

K=Data['K']
M=Data['M']

Theta=np.arctan(M/K) # phase speed angle. dispersion is omega^2=N^2 cos(Theta)^2 

theta=np.arange(51)*np.pi/100 # ranges of [0,pi/2]

Usp=Data['Usp']
time=Data['time']
dt=np.mean(time[1:]-time[:-1])
window=np.hamming(128)[:,np.newaxis,np.newaxis]/128
nr=int(np.floor(time.shape[0]/128))-1
print("calculating int E(omega,k cos(theta),k sin(theta))kdk")
SpReduced=np.zeros((2*nr,64,51))
with Bar("Processing ", max=2*nr) as bar:
    for n in range(2*nr):
        data=Usp[n*64:n*64+128,:,:]*window
        sp=fft(data,axis=0)
        SP=np.abs(sp[:64,:,:])**2
        for k in range(51):
            for o in range(64):
                dummy=SP[o,:,:]
                SpReduced[n,o,k]=np.sum(dummy[Theta>=theta[k]])
        bar.next()
SpReduced=SpReduced[:,:,:-1]-SpReduced[:,:,1:]

SPOM=np.zeros((2*nr,64, Usp.shape[2]))
print("calculating int E(omega,k,m) dk")
with Bar("Processing ", max=2*nr) as bar:
    for n in range(2*nr):
        data=Usp[n*64:n*64+128,:,:]*window
        sp=fft(data,axis=0)
        SP=np.abs(sp[:64,:,:])**2
        SPOM[n,:,:]=np.sum(SP,axis=1)
        bar.next()

SPOK=np.zeros((2*nr,64, Usp.shape[1]))
print("calculating int E(omega,k,m) dm")
with Bar("Processing ", max=2*nr) as bar:
    for n in range(2*nr):
        data=Usp[n*64:n*64+128,:,:]*window
        sp=fft(data,axis=0)
        SP=np.abs(sp[:64,:,:])**2
        SPOK[n,:,:]=np.sum(SP,axis=2)
        bar.next()







omega=np.arange(64)*2*np.pi/(dt*128) # we assume that Dt between saves is 2. 
np.savez(outputfile, SP=SpReduced, SPOM=SPOM, SPOK=SPOK,Omega=omega, theta=(theta[1:]+theta[:-1])/2)


