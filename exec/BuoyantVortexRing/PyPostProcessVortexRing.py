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
sys.path.append(os.environ['HOME']+'/LES-SOMAR/PythonScripts/')
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



def calc_CM(inputDir, outputDir, beg=0, end=None, step=None, level=0):
    import os
    from tqdm import tqdm
    import matplotlib.pyplot as plt
    import numpy as np
    os.chdir(inputDir)
    
    print(outputDir)
    CWD = os.getcwd()
    filenames = [f for f in os.listdir(CWD) if os.path.isfile(os.path.join(CWD, f)) and f[0:4] == 'plot' and f[-8:] == '.3d.hdf5']
    filenames.sort()
    filenames=filenames[beg:end:step]
    try:
        Nx=GetNx(filenames[0],level)
    except:
        Nx=GetNx(filenames[0],level-1)
    CM=np.zeros((len(filenames),4))
    
    time=np.zeros(len(filenames))
    
    parmDict=getParms(inputname='/home/adscotti/LES-SOMAR/exec/VortexRing/inputs')
    
    Lx,Ly, Lz=parmDict['base.L']
    
    n=0
    X=np.zeros((len(filenames),))
    Y=np.zeros((len(filenames),))
    Z=np.zeros((len(filenames),))

    for filename in tqdm(filenames):
        try:
            Q, dx, time[n] = GatherComponents(filename, comps=('T_total',),level=level)
        except:
            Q, dx, time[n] = GatherComponents(filename, comps=('T_total',),level=level-1)
        b=-Q[:,:,:,0]
        bm=np.nanmean(b)
        r2 = ((np.arange(-Nx[0]//2,Nx[0]//2)[:,np.newaxis, np.newaxis]+.5) * dx[0])**2 + ((np.arange(-Nx[1]//2,Nx[1]//2)[np.newaxis:,np.newaxis]+.5) * dx[1])**2 
        CM[n,2] = np.nanmean(b * (np.arange(-Nx[2]//2,Nx[2]//2)[np.newaxis,np.newaxis,:]+.5) * dx[2])/bm
        CM[n,3] = np.sqrt(np.nanmean(b * r2)/bm) 
        bm = np.nanmean(np.squeeze(b[:,Nx[1]//2-1,:]))
        CM[n,0] = np.nanmean(np.squeeze(b[:,Nx[1]//2-1,:]) * np.abs((np.arange(-Nx[0]//2,Nx[0]//2)[:,np.newaxis]+.5) * dx[0]))/(2*bm)
        bm = np.nanmean(b[Nx[0]//2-1,:,:])
        CM[n,1] = np.nanmean(np.squeeze(b[Nx[0]//2-1,:,:]) * np.abs((np.arange(-Nx[1]//2,Nx[1]//2)[:,np.newaxis]+.5) * dx[1]))/(2*bm)
        
        
        n+=1


    plt.plot( time, CM[:,0], label='x')
    plt.plot( time, CM[:,1], label='y')
    plt.plot( time, CM[:,2], label='z')
    plt.plot(time, CM[:,3], label='R')
    plt.xlabel('Time')
    plt.ylabel('CM')
    plt.title('Center of Mass vs Time')
    plt.legend()
    plt.grid()
    plt.show()
    
    np.savez(outputDir+'/results',time=time, CM=CM)




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
    beg = 0
    end = None
    step = 1
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
    
    
    calc_CM(inputDir, outputDir, beg=beg, end=end, step=step,level=2)
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










