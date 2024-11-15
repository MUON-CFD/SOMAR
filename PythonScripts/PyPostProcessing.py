#!/usr/bin/env python3
import re
import numpy as np
import sys
from scipy.fftpack import dct, dst 
sys.path.append('/home/alberto/SOMAR_V1.1/LES-SOMAR/PythonScripts/')

import ChomboIO as ChIO

from Box import Box
import matplotlib.pyplot as plt

def GetNx(filename, level=0):
    Fh=ChIO.ReadCheckPoint(filename)
    if level>Fh.RootAttributes['max_level']:
        raise ValueError("Level requested does not exist")
    domain=Fh.LevelsAttributes[level]['prob_domain']
    NX=[]
    dim=len(domain)//2
    for d in range(dim):
        NX.append(domain[dim+d]-domain[d]+1)
    del Fh
    return NX

def GatherComponents(filename, comps, level=0):
    Fh = ChIO.ReadCheckPoint(filename)
    if level > Fh.RootAttributes['max_level']:
        raise ValueError("Level requested does not exist")

    LDFAB = Fh.SetUpLevelDatasFABs(FABName='data')
    Fh.ReadFABs(LevelDatas=LDFAB, FABName='data')

    box = tuple(Fh.LevelsAttributes[level]['prob_domain'])

    b= Box(box+('Box',))

    CompNumber=list()
    for name in comps:

        for (key, value) in Fh.RootAttributes.items():
             if type(value) == np.bytes_:

                 if str(value)[2:-1] == name:
                    CompNumber.append(int(key[10:]))


    if len(CompNumber) != len(comps):
        raise ValueError("PRoblem with comps names")

    X = LDFAB[level].GatherBox(b, components=CompNumber)

    return X, Fh.LevelsAttributes[level]['vec_dx'], Fh.RootAttributes['time']

def contourPlot(x,y, zrange,Q,comp=0,imagetitle='title',filename="test.png"):
    xm,xM=x
    ym,yM=y
    nx,ny=Q[:,:,comp].shape
    dx=(xM-xm)/nx
    dy=(yM-ym)/ny
    z_min,z_max=zrange

    #X,Y=np.mgrid(slice(-ym,yM+dy,dy), slice(-xm,xM+dx,dx))
    plt.imshow(Q[:,:,comp].transpose(), cmap='RdBu', vmin=z_min, vmax=z_max,
           extent=[xm, xM, ym, yM])
    plt.title(imagetitle)
    plt.colorbar()
    plt.savefig(filename, format='png')
    plt.close()

def getParms(inputname):
    parmDict={}
    f=open(inputname,'r')
    parms=f.read()
    f.close()
    ss=parms.split('\n') # split lines
    x=[s.split() for s in ss if (len(s)>0) and s[0]!='#'] # eliminate comment lines
    for line in x:
        y=[]
        j=[i for i,entry in enumerate(line) if entry=='#'] # only consider entries up to a #
        if len(j)>0:
            ll=line[2:j[0]]
        else:
            ll=line[2:]
        for entry in ll:
            try:
                y.append(int(entry))
            except:
                try:
                    y.append(float(entry))
                except:
                    y.append(entry)
            finally:
                parmDict[line[0]] = y[0] if len(y)==1 else tuple(y)

    return parmDict

if __name__ == "__main__":

    import os

    os.chdir('../execGarrettMunk/hdf5_output/')
    from progress.bar import Bar

    CWD = os.getcwd()
    filenames = [f for f in os.listdir(CWD) if os.path.isfile(os.path.join(CWD, f)) and f[0:4] == 'plot']
    filenames.sort()
    Nx=GetNx(filenames[0])
    BF=np.zeros(len(filenames))
    prod = np.zeros(len(filenames))
    dissAPE=np.zeros(len(filenames))
    diss = np.zeros(len(filenames))
    KE = np.zeros(len(filenames))
    APE=np.zeros(len(filenames))
    KEFlux = np.zeros(len(filenames))
    RiRT = np.zeros(len(filenames))
    RiKH=np.zeros_like(RiRT)
    time = np.zeros(len(filenames))
    U = np.zeros((len(filenames),Nx[1]))
    V = np.zeros((len(filenames),Nx[1]))
    B = np.zeros((len(filenames), Nx[1]))
    BAve=np.zeros((len(filenames), Nx[1]))
    Usp=np.zeros((len(filenames), Nx[0], Nx[1]))
    K=np.zeros(Nx)
    M=np.zeros(Nx)
    omega=np.zeros(Nx)
    parmDict=getParms(inputname='../inputs.template.machine')
    #for key,value in parmDict.items():
    #    print(key," : ", value)
    nu=parmDict['rhs.nu']
    kappa=parmDict['rhs.bKappa']
    Lx,Ly=parmDict['rhs.L']
    Omega=np.sin(parmDict['rhs.theta'])
    Nsq=parmDict['rhs.Nsq']
    K[:,:]=np.pi/Lx*(np.arange(Nx[0])+1)[:,np.newaxis]*np.ones(Nx[1])[np.newaxis,:]
    M[:,:]=np.ones(Nx[0])[:,np.newaxis]*np.pi/Ly*np.arange(Nx[1])[np.newaxis,:]
    omega=np.sqrt(Nsq*K**2/(K**2+M**2))
    with Bar('Processing', max=len(filenames)) as bar:
        for n,filename in enumerate(filenames):


            Q, dx, time[n] = GatherComponents(filename, comps=('x_vel', 'y_vel', 'b_pert', 'pressure'))
            u=Q[:,:,0]
            Usp[n,:,:]=dst(dct(u,type=2,axis=1),type=2,axis=0)

            timeOmega = time[n]*Omega / (2 * np.pi)
            contourPlot((0,30),(-5,5),(-.1,.1),Q,comp=2,imagetitle=str(timeOmega),filename=str(n)+'.png')
            ke=(Q[:, :, 0]** 2 + Q[:, :, 1]** 2) / 2
            KE[n] = np.sum(ke) * dx[0] * dx[1]
            ape=(Q[:,:,2]**2)/2
            APE[n] = np.sum(ape) *dx[0] *dx[1]
            BF[n] = np.sum((Q[:,:,1]*Q[:,:,2])) *dx[0] *dx[1] # note this does not include contribution
                                                             # from w*bbar. The latter is zero only averaged over cycles.
            S11 = ((Q[2:   , 1: -1, 0] - Q[0: -2, 1: -1, 0]) / (2 * dx[0]))
            S22 = ((Q[1: -1, 2:   , 1] - Q[1: -1, 0: -2, 1]) / (2 * dx[1]))
            S12 = .5 * ((Q[1:-1, 2: , 0] - Q[1:-1, 0:-2, 0]) / (2 * dx[1]) +
                        (Q[2: , 1:-1, 1] - Q[0:-2, 1:-1, 1]) / (2 * dx[0]))

            G1 = ((Q[2:, 1:-1, 2] - Q[0:-2, 1:-1, 2]) / (2 * dx[0]))
            G2 = ((Q[1:-1, 2:, 2] - Q[1:-1, 0:-2, 2]) / (2 * dx[1]))
            epsilonAPE=nu*(G1**2+G2**2)
            epsilon = nu * (2 * S11 ** 2 + 2 * S22 ** 2 + 4*S12 ** 2)
            pu = Q[0, :, 0] * Q[0, :, 3]
            KEFlux[n] = np.sum(Q[0,:,0]*ke[0,:]) * dx[1]
            prod[n] = np.sum(pu) *dx[1]
            diss[n] = np.sum(epsilon) * dx[0] * dx[1]
            dissAPE[n] = np.sum(epsilonAPE) * dx[0] * dx[1]

            S2 = ((Q[:, 2:, 0] - Q[:, 0:-2,0]) / (2*dx[1]))** 2

            N2 =1.0 -(Q[:, 2:, 2] - Q[:, 0:-2,2]) / (2 * dx[1])

            Ri = N2 / S2

            Count, Edges = np.histogram(Ri, bins=[min(Ri.min(), -1),0, .25, 1, Ri.max()])

            RiRT[n] = Count[0] / (Count[0] + Count[1] + Count[2] + Count[3])
            RiKH[n] = Count[1] / (Count[0] + Count[1] + Count[2] + Count[3])
            U[n,:] = Q[128,:, 0]
            V[n,:] = Q[128,:, 1]
            B[n,:] = Q[128,:,2]
            BAve[n,:]=np.sum(Q[:,:,2],axis=0)/Q.data.shape[0]
            bar.next()


    np.savez('results',time=time, K=K, M=M, omega=omega,Usp=Usp, KE=KE, APE=APE, RiInst=RiRT, RiKH=RiKH, U=U, V=V, B=B, BAve=BAve, diss=diss, prod=prod, dissAPE=dissAPE, BF=BF, KEFlux=KEFlux)
    framerate=24
    try:
        os.system("ffmpeg -r " +str(framerate) +" -i %01d.png -vcodec mpeg4 -y movie.mp4")
        os.system("rm *.png")
    except:
        pass







