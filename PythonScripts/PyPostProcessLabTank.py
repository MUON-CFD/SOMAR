#!/usr/bin/env python3
from lib2to3.pgen2.pgen import generate_grammar
from logging import root
import re
import numpy as np
import sys
from scipy.fftpack import dct, dst 
from scipy.integrate import cumtrapz
from scipy.io import savemat


import os
sys.path.append(os.environ['HOME']+'/SOMAR_V1.1/LES-SOMAR/PythonScripts/')
from mpi4py import MPI
import ChomboIO as ChIO

from Box import Box
import matplotlib.pyplot as plt
import Diagnostic 
from mpi4py import MPI
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

def zstar(inFilename, outFilename):
    from LevelDataFAB import LevelDataFAB as LDFAB
    Fh = ChIO.ReadCheckPoint(inFilename)
    

    dataLDFABs = Fh.SetUpLevelDatasFABs(FABName='data')
    Fh.ReadFABs(LevelDatas=dataLDFABs, FABName='data')
    Dx=Fh.LevelsAttributes[0]['vec_dx']
    L=Fh.LevelsAttributes[0]['prob_domain']
    SpaceDim=len(Dx)
    dV=1
    
    
    for dx in Dx:
        dV*=dx

    A=1
    for dim in range(SpaceDim-1):
    
        A*=(L[SpaceDim+dim]-L[dim])*Dx[dim]
        
    CompNumber=list()
    for name in ['b_total']:
        for (key, value) in Fh.RootAttributes.items():
             if type(value) == np.bytes_:

                 if str(value)[2:-1] == name:
                    CompNumber.append(int(key[10:]))

    hasBTotal = bool(len(CompNumber)==1)
    
    if not hasBTotal:
        for name in ['S_pert']:
            for (key, value) in Fh.RootAttributes.items():
                if type(value) == np.bytes_:

                    if str(value)[2:-1] == name:
                        CompNumber.append(int(key[10:]))
        import h5py as h5
        F=h5.File('plot_stratData_lev0.hdf5')
        nz=F['level_0'].attrs['prob_domain'][2*SpaceDim-1]-F['level_0'].attrs['prob_domain'][SpaceDim-1] + 1
        zOffset=F['level_0'].attrs['prob_domain'][SpaceDim-1]
        dataOffset=(nz+2) * (3 ** (SpaceDim-1))
        k = 1 # bbar is store in component 1
        if SpaceDim == 2:
            bbar=F['level_0/data:datatype=0'][k * dataOffset: (k+1) * dataOffset].reshape(nz+2,3)[1:-1,1]
        else:
            bbar=F['level_0/data:datatype=0'][k * dataOffset: (k+1) * dataOffset].reshape(nz+2,3,3)[1:-1,1,1]



    LD=dataLDFABs[0]
    M=-1e6
    m=1e6
    
    for _data in LD:
        if hasBTotal:
            try:
                M=max(M,np.amax(_data[1:-1,1:-1,1:-1,CompNumber[0]]))
                m=min(m,np.amin(_data[1:-1,1:-1,1:-1,CompNumber[0]]))
            except:
                M=max(M,np.amax(_data[1:-1,1:-1,CompNumber[0]]))
                m=min(m,np.amin(_data[1:-1,1:-1,CompNumber[0]]))

        else:
            try:
                bTotal=np.empty_like(_data[1:-1,1:-1,1:-1,CompNumber[0]])
                bTotal[...]=_data[1:-1,1:-1,1:-1,CompNumber[0]]
            except:
                bTotal=np.empty_like(_data[1:-1,1:-1,CompNumber[0]])
                bTotal[...]=_data[1:-1,1:-1,CompNumber[0]]

            box=Box.FromAnotherBox(_data.box)
            
            box.Grow(SpaceDim*[-1])
            #Z=np.empty(box.Size())
            
            #Z[...]=np.arange(box.LoEnd()[-1],box.HiEnd()[-1]+1)*Dx[-1]+Dx[-1]/2

            #bTotal[...]+=-355*(Z[...]-0.43)
            bTotal[...]*=9.81*8e-4
            bTotal[...]+=-9.81*2e-4*20

            
            bTotal[...]+=bbar[box.LoEnd()[SpaceDim-1]-zOffset:box.HiEnd()[SpaceDim-1]-zOffset+1]



            M=max(M,np.amax(bTotal))
            m=min(m,np.amin(bTotal))




            

    MPI.COMM_WORLD.Barrier()
    MPI.COMM_WORLD.allreduce(M,op=MPI.MAX)
    MPI.COMM_WORLD.allreduce(m,op=MPI.MIN)

    
    N=np.zeros(200)
    for _data in LD:
        if hasBTotal:
            try:
                n,edges = np.histogram(_data[1:-1,1:-1,1:-1,CompNumber[0]], bins=200, range=(m,M))
            except:
                n,edges = np.histogram(_data[1:-1,1:-1,CompNumber[0]], bins=200, range=(m,M))
        else:
            try:
                bTotal=np.empty_like(_data[1:-1,1:-1,1:-1,CompNumber[0]])
                bTotal[...]=_data[1:-1,1:-1,1:-1,CompNumber[0]]
            except:
                bTotal=np.empty_like(_data[1:-1,1:-1,CompNumber[0]])
                bTotal[...]=_data[1:-1,1:-1,CompNumber[0]]

            box=Box.FromAnotherBox(_data.box)
            
            box.Grow(SpaceDim*[-1])
            #Z=np.empty(box.Size())
            
            #Z[...]=np.arange(box.LoEnd()[-1],box.HiEnd()[-1]+1)*Dx[-1]+Dx[-1]/2

            #bTotal[...]+=-355*(Z[...]-0.43)
            bTotal[...]*=9.81*8e-4
            bTotal[...]+=-9.81*2e-4*20
            bTotal[...]+=bbar[box.LoEnd()[SpaceDim-1]-zOffset:box.HiEnd()[SpaceDim-1]-zOffset+1]
            n,edges = np.histogram(bTotal, bins=200, range=(m,M))

        N+=n*dV
    
    if MPI.COMM_WORLD.rank==0:
        Ntot=np.zeros_like(N)
    else:
        Ntot=None

    MPI.COMM_WORLD.Barrier()
    MPI.COMM_WORLD.Reduce(N,Ntot,op=MPI.SUM,root=0)
    MPI.COMM_WORLD.Barrier()

    if MPI.COMM_WORLD.rank==0:
        zs=np.cumsum(Ntot)/A

        data=dict()
        data['zs']=zs
        data['bs']=(edges[:-1]+edges[1:])/2
        savemat(outFilename,data)
    return 

def GammaEfficient(infile, N2, nu=1e-6):
    from LevelDataFAB import LevelDataFAB as LDFAB
    Fh=ChIO.ReadCheckPoint(infile)
    qLDFABs = Fh.SetUpLevelDatasFABs(FABName='data')
    Fh.ReadFABs(LevelDatas=qLDFABs, FABName='data')
    gamma = 0
    chiTotal=0
    epsilonTotal=0
    TKETotal = 0
    time = Fh.RootAttributes['time']
    Comp=dict()
    for name in ['epsilon', 'chi', 'x_vel', 'y_vel', 'z_vel']:
        
        for (key, value) in Fh.RootAttributes.items():
             
             if type(value) == np.bytes_:

                 if str(value)[2:-1] == name:
                     Comp[name]=int(key[10:])

    
    if len(qLDFABs)>1: 
        qLDFAB = qLDFABs[1]
        
    else:
        # chiLDFAB = chiLDFABs[0]
        # epsilonLDFAB = epsilonLDFABs[0]

        return gamma, 0, 0, 0, time
    for _q in qLDFAB:
        chiTotal += _q[...,Comp['chi']].mean()
        epsilonTotal += _q[...,Comp['epsilon']].mean()
        TKETotal += ((_q[...,Comp['x_vel']]**2).mean() +
                     (_q[...,Comp['y_vel']]**2).mean() +
                      (_q[...,Comp['z_vel']]**2).mean())/2

    APETotal = chiTotal/(2*N2)
    
    gamma = APETotal/(APETotal+epsilonTotal)


    Rib = epsilonTotal/(nu * N2)
    
    Fr = epsilonTotal/(np.sqrt(N2) * TKETotal)

    return gamma, Rib, Fr, epsilonTotal, time


def Gamma(inChi, inDiss, inQ, N2, nu = 1e-6):
    from LevelDataFAB import LevelDataFAB as LDFAB
    FhDiss = ChIO.ReadCheckPoint(inDiss)
    FhChi = ChIO.ReadCheckPoint(inChi)
    FhQ = ChIO.ReadCheckPoint(inQ)

    epsilonLDFABs = FhDiss.SetUpLevelDatasFABs(FABName='data')
    FhDiss.ReadFABs(LevelDatas=epsilonLDFABs, FABName='data')

    chiLDFABs = FhChi.SetUpLevelDatasFABs(FABName='data')
    FhChi.ReadFABs(LevelDatas=chiLDFABs, FABName='data')

    qLDFABs = FhQ.SetUpLevelDatasFABs(FABName='data')
    FhQ.ReadFABs(LevelDatas=qLDFABs, FABName='data')

    gamma = 0
    chiTotal=0
    epsilonTotal=0
    TKETotal = 0
    time = FhDiss.RootAttributes['time']
    if len(chiLDFABs)>1: 
        chiLDFAB = chiLDFABs[1]
        epsilonLDFAB = epsilonLDFABs[1]
        qLDFAB = qLDFABs[1]
    else:
        # chiLDFAB = chiLDFABs[0]
        # epsilonLDFAB = epsilonLDFABs[0]
        return gamma, 0, 0, time
        
    for (_c,_d,_q)  in zip(chiLDFAB, epsilonLDFAB,qLDFAB):
        
        # we consider interior points to avoid issues with CF boundaries
        if _d.dimension==3:
            chiTotal +=_c.data[1:-1,1:-1,1:-1].mean()
            epsilonTotal +=_d.data[1:-1,1:-1,1:-1].mean()
            TKETotal += 3/2 * (_q.data[1:-1,1:-1,1:-1,1].mean())**2
        else:
            chiTotal +=_c.data[1:-1,1:-1].mean()
            epsilonTotal +=_d.data[1:-1,1:-1].mean()
            TKETotal += 3/2 * (_q.data[1:-1,1:-1,1].mean())**2
    
    APETotal = chiTotal/(2*N2)
    
    gamma = APETotal/(APETotal+epsilonTotal)


    Rib = epsilonTotal/(nu * N2)
    
    Fr = epsilonTotal/(np.sqrt(N2) * TKETotal)

    return gamma, Rib, Fr, time

def Chi(k,inFilename, outFilename):
    from LevelDataFAB import LevelDataFAB as LDFAB
    Fh = ChIO.ReadCheckPoint(inFilename)
    

    dataLDFABs = Fh.SetUpLevelDatasFABs(FABName='data')
    Fh.ReadFABs(LevelDatas=dataLDFABs, FABName='data')

    CompNumber=list()
    for name in ['eddyNu']:
        
        for (key, value) in Fh.RootAttributes.items():
             if type(value) == np.bytes_:

                 if str(value)[2:-1] == name:
                    CompNumber.append(int(key[10:]))

    for name in ['S_pert']:
        
        for (key, value) in Fh.RootAttributes.items():
             if type(value) == np.bytes_:

                 if str(value)[2:-1] == name:
                    CompNumber.append(int(key[10:]))

    ChiLDFABs = tuple([LDFAB.MakeFromLDFAB(dataLDFAB) for dataLDFAB in dataLDFABs ])
    for level,(_eLDFAB,_dataLDFAB) in enumerate(zip(ChiLDFABs,dataLDFABs)):
        _eLDFAB.nComp = 1
        for (_e,_d) in zip(_eLDFAB,_dataLDFAB):
            DX=tuple(Fh.LevelsAttributes[level]['vec_dx'])+('RealVect',)
            
            _e.data=np.empty_like(Diagnostic.differentiateFAB(_d,DX,CompNumber[1],0))
            for d in range(_d.dimension):
                dummy=9.8*8e-4*Diagnostic.differentiateFAB(_d,DX,CompNumber[1],d)
                
                
                try:
                    if d>0:
                        _e.data+=dummy*dummy*(k+_d[1:-1,1:-1,1:-1,CompNumber[0]])
                    else:
                        _e.data[...]=dummy*dummy*(k+_d[1:-1,1:-1,1:-1,CompNumber[0]])
                except:
                    if d>0:
                        _e.data+=dummy*dummy*(k+_d[1:-1,1:-1,CompNumber[0]])
                    else:
                        _e.data[...]=dummy*dummy*(k+_d[1:-1,1:-1,CompNumber[0]])
                
            _e.data.reshape((_e.data.shape)+(1,))
            _e.box.Grow((-1,-1,-1))
            

    W=ChIO.WriteCheckPoint(outFilename)

    W.WriteHeader('Chombo_global','SpaceDim',dataLDFABs[0].SpaceDim)
    W.WriteHeader('Chombo_global', 'testReal',float(0))

    W.WriteHeader('/', 'num_components', 1)
    W.WriteHeader('/', 'num_levels', Fh.RootAttributes['num_levels'])
    W.WriteHeader('/', 'time', Fh.RootAttributes['time'])

    W.WriteHeader('/' , 'max_level', Fh.RootAttributes['max_level'])

    W.WriteHeader('/', 'component_0', np.string_('Chi'))
    W.WriteHeader('/', 'iteration', Fh.RootAttributes['iteration'])

    for level,(_eLDFAB,attrs) in enumerate(zip(ChiLDFABs,Fh.LevelsAttributes)):
        W.WriteHeader('level_'+str(level), 'dt', attrs['dt'])
        W.WriteHeader('level_'+str(level), 'prob_domain', attrs['prob_domain'])
        W.WriteHeader('level_'+str(level), 'time', attrs['time'])
        W.WriteHeader('level_'+str(level), 'ref_ratio', attrs['ref_ratio'])
        W.WriteHeader('level_'+str(level), 'vec_dx', attrs['vec_dx'])
        W.WriteLevelDataFAB(_eLDFAB, 'level_'+str(level), 'data')


    del W
    del Fh

def dissipation(nu,inFilename, outFilename):
    from LevelDataFAB import LevelDataFAB as LDFAB
    Fh = ChIO.ReadCheckPoint(inFilename)
    

    dataLDFABs = Fh.SetUpLevelDatasFABs(FABName='data')
    Fh.ReadFABs(LevelDatas=dataLDFABs, FABName='data')

    CompNumber=list()
    for name in ['eddyNu']:
        
        for (key, value) in Fh.RootAttributes.items():
             if type(value) == np.bytes_:

                 if str(value)[2:-1] == name:
                    CompNumber.append(int(key[10:]))
    
    epsilonLDFABs = tuple([LDFAB.MakeFromLDFAB(dataLDFAB) for dataLDFAB in dataLDFABs ])

    for level,(_eLDFAB,_dataLDFAB) in enumerate(zip(epsilonLDFABs,dataLDFABs)):
        _eLDFAB.nComp = 1
        for (_e,_d) in zip(_eLDFAB,_dataLDFAB):
            DX=tuple(Fh.LevelsAttributes[level]['vec_dx'])+('RealVect',)
            
            S=Diagnostic.Strain(_d,DX)
             
                
            try:
                nut=(nu + _d[1:-1,1:-1,1:-1,CompNumber[0]])
            except:
                nut=(nu + _d[1:-1,1:-1,CompNumber[0]])
            _e.data=Diagnostic.traceStrainSquared(S)* 2 * nut 
            
            _e.data.reshape((_e.data.shape)+(1,))
            _e.box.Grow((-1,-1,-1))
            

    W=ChIO.WriteCheckPoint(outFilename)

    W.WriteHeader('Chombo_global','SpaceDim',dataLDFABs[0].SpaceDim)
    W.WriteHeader('Chombo_global', 'testReal',float(0))

    W.WriteHeader('/', 'num_components', 1)
    W.WriteHeader('/', 'num_levels', Fh.RootAttributes['num_levels'])
    W.WriteHeader('/', 'time', Fh.RootAttributes['time'])

    W.WriteHeader('/' , 'max_level', Fh.RootAttributes['max_level'])

    W.WriteHeader('/', 'component_0', np.string_('diss'))
    W.WriteHeader('/', 'iteration', Fh.RootAttributes['iteration'])

    for level,(_eLDFAB,attrs) in enumerate(zip(epsilonLDFABs,Fh.LevelsAttributes)):
        W.WriteHeader('level_'+str(level), 'dt', attrs['dt'])
        W.WriteHeader('level_'+str(level), 'prob_domain', attrs['prob_domain'])
        W.WriteHeader('level_'+str(level), 'time', attrs['time'])
        W.WriteHeader('level_'+str(level), 'ref_ratio', attrs['ref_ratio'])
        W.WriteHeader('level_'+str(level), 'vec_dx', attrs['vec_dx'])
        W.WriteLevelDataFAB(_eLDFAB, 'level_'+str(level), 'data')


    del W
    del Fh








    


    

    


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

def calc_spectra(inputDir, outputDir, beg=0, end=None, step=None):
    import os

    os.chdir(inputDir)
    from progress.bar import Bar
    print(outputDir)
    CWD = os.getcwd()
    filenames = [f for f in os.listdir(CWD) if os.path.isfile(os.path.join(CWD, f)) and f[0:4] == 'plot' and f[-8:] == '.2d.hdf5']
    filenames.sort()
    filenames=filenames[beg:end:step]
    Nx=GetNx(filenames[0])
    Usp=np.zeros((len(filenames), Nx[0], Nx[1]))
    K=np.zeros(Nx)
    M=np.zeros(Nx)
    B=np.zeros(len(filenames))

    time=np.zeros(len(filenames))
    omega=np.zeros(Nx)
    parmDict=getParms(inputname='/home/adscotti/LES-SOMAR/exec/LabTank/inputs2D')
    #for key,value in parmDict.items():
    #    print(key," : ", value)
    
    Lx,Ly=parmDict['rhs.L']
    
    Nsq=parmDict['rhs.dSdZbar']*9.81*(-8e-4) 
    Omega=parmDict['forcing.frequency'] 
    K[:,:]=np.pi/Lx*(np.arange(Nx[0])+1)[:,np.newaxis]*np.ones(Nx[1])[np.newaxis,:]
    M[:,:]=np.ones(Nx[0])[:,np.newaxis]*np.pi/Ly*np.arange(Nx[1])[np.newaxis,:]
    omega=np.sqrt(Nsq*K**2/(K**2+M**2))
    with Bar('Processing', max=len(filenames)) as bar:
        for n,filename in enumerate(filenames):
            Q, dx, time[n] = GatherComponents(filename, comps=('x_vel', 'S_pert'))
            
            u=Q[:,:,0]
            Usp[n,:,:]=dst(dct(u,type=2,axis=1),type=2,axis=0)
            B[n]=Q[530,1800,1]*(-8e-4*9.81)
            timeOmega = time[n]*Omega / (2 * np.pi)
            bar.next()
    np.savez(outputDir+'/results',time=time, K=K, M=M, omega=omega,Usp=Usp, b_pert=B)
    return

def main():
    import os

    os.chdir('../exec/LabTank/hdf5_output/')
    from progress.bar import Bar

    CWD = os.getcwd()
    filenames = [f for f in os.listdir(CWD) if os.path.isfile(os.path.join(CWD, f)) and f[0:4] == 'plot' and f[-8:] == '.2d.hdf5']
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
    parmDict=getParms(inputname='../inputs2D')
    #for key,value in parmDict.items():
    #    print(key," : ", value)
    nu=parmDict['rhs.nu']
    kappa=parmDict['rhs.SKappa']
    Lx,Ly=parmDict['rhs.L']
    
    Nsq=parmDict['rhs.dSdZbar']*9.81*(-8e-4) 
    Omega=np.sin(parmDict['forcing.theta']) *np.sqrt(Nsq)
    K[:,:]=np.pi/Lx*(np.arange(Nx[0])+1)[:,np.newaxis]*np.ones(Nx[1])[np.newaxis,:]
    M[:,:]=np.ones(Nx[0])[:,np.newaxis]*np.pi/Ly*np.arange(Nx[1])[np.newaxis,:]
    omega=np.sqrt(Nsq*K**2/(K**2+M**2))
    with Bar('Processing', max=len(filenames)) as bar:
        for n,filename in enumerate(filenames):


            Q, dx, time[n] = GatherComponents(filename, comps=('x_vel', 'y_vel', 'S_pert', 'pressure'))
            Q[:,:,2] *= 8e-4*9.81 # b_pert
            if n>0:
                print(" ", time[n]-time[n-1])
            
            u=Q[:,:,0]
            Usp[n,:,:]=dst(dct(u,type=2,axis=1),type=2,axis=0)

            timeOmega = time[n]*Omega / (2 * np.pi)
            contourPlot((0,30),(-5,5),(-.1,.1),Q,comp=2,imagetitle=str(timeOmega),filename=str(n)+'.png')
            ke=(Q[:, :, 0]** 2 + Q[:, :, 1]** 2) / 2
            KE[n] = np.sum(ke) * dx[0] * dx[1]
            ape=(Q[:,:,2]**2)/2/Nsq
            APE[n] = np.sum(ape) *dx[0] *dx[1]
            BF[n] = np.sum((Q[:,:,1]*Q[:,:,2])) *dx[0] *dx[1] # note this does not include contribution
                                                             # from w*bbar. The latter is zero only averaged over cycles.
            S11 = ((Q[2:   , 1: -1, 0] - Q[0: -2, 1: -1, 0]) / (2 * dx[0]))
            S22 = ((Q[1: -1, 2:   , 1] - Q[1: -1, 0: -2, 1]) / (2 * dx[1]))
            S12 = .5 * ((Q[1:-1, 2: , 0] - Q[1:-1, 0:-2, 0]) / (2 * dx[1]) +
                        (Q[2: , 1:-1, 1] - Q[0:-2, 1:-1, 1]) / (2 * dx[0]))

            G1 = ((Q[2:, 1:-1, 2] - Q[0:-2, 1:-1, 2]) / (2 * dx[0]))
            G2 = ((Q[1:-1, 2:, 2] - Q[1:-1, 0:-2, 2]) / (2 * dx[1]))
            epsilonAPE=kappa*(G1**2+G2**2)
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

if __name__ == "__main__":
    import getopt
    options="i:o:b:e:s:h"
    long_options=['help', 'inputDir', 'outputDir', 'begin', 'end', 'step']
    inputDir=''
    outputDir=''
    try:
        arguments, values = getopt.getopt(sys.argv[1:],options,long_options)
        
    except getopt.GetoptError:
        if MPI.COMM_WORLD.Get_rank()==0:
            print("Usage: PyPostProcessing.py -i <inputDir> -o <outputDir> -b first sequential file -e last sequential file ")
        MPI.COMM_WORLD.Barrier()
        sys.exit(2)
    
    for argument, value in arguments:
        if argument == '-i':
            inputDir=value
        if argument == '-o':
            outputDir = value
        if argument == '-b':
            beg=int(value)
        if argument == '-e':
            end=int(value)
        if argument == '-s':
            step=int(value)
        if argument == '-h':
            if MPI.COMM_WORLD.Get_rank()==0:
                print("Usage: mpirun -np XX PyPostProcessing.py -i <inputDir> -o <outputDir> \n")
                print("mpirun -np xx is necessary only if multiple cores are requested")
            MPI.COMM_WORLD.Barrier()
            sys.exit(0)
    
    calc_spectra(inputDir, outputDir,beg=beg,end=end,step=step)
    
    # os.chdir(inputDir)
    # from progress.bar import Bar

    # CWD = os.getcwd()
    # filenames = [f for f in os.listdir(CWD) if os.path.isfile(os.path.join(CWD, f)) and f[0:4] == 'plot' and (f[-8:] == '.2d.hdf5' or f[-8:] == '.3d.hdf5')]
    # filenames.sort()
    # gamma=[]
    # Reb = []
    # Fr = []
    # time = []
    # epsilon = []
    # with Bar('Processing', max=len(filenames)) as bar:
    #     for n,filename in enumerate(filenames):
    #         #dissipation(1e-6, filename, 'diss'+filename[4:])
    #         #zstar(filename, 'zs'+filename[4:11]+'.mat' )
    #         #Chi(1e-7, filename, 'chi'+filename[4:])
    #         #try:
    #         #Rif,I,fr,t  = Gamma('chi'+filename[4:],'diss'+filename[4:], filename, 9.81*8e-4*335)
    #         Rif,I,fr, e, t = GammaEfficient(filename, 9.81 * 8e-4 * 16)
    #         gamma.append(Rif)
    #         Reb.append(I)
    #         Fr.append(fr)
    #         epsilon.append(e)
    #         time.append(t)
    #         #except:
    #         #    pass
    #         #if n>500: break
    #         if MPI.COMM_WORLD.rank==0:            
    #             bar.next()

    # dummy={ 'gamma' : gamma, 'Reb' : Reb, 'Fr' : Fr, 'epsilon' : epsilon, 'time' : time}

    # savemat('gamma.mat', dummy)










