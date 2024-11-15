#!/usr/bin/env python3
SOMARSCRIPTDIR='/home/alberto/SOMAR_V1.1/LES-SOMAR/PythonScripts'
import re
import numpy as np
from scipy import integrate as integrate
import sys
from scipy.fftpack import dct, dst
sys.path.append(SOMARSCRIPTDIR)
from mpi4py import MPI
import ChomboIO as ChIO

from Box import Box
import matplotlib.pyplot as plt
import Diagnostic
from LevelDataFAB import LevelDataFAB
def GatherComponents(filename, level=0):
    Fh = ChIO.ReadCheckPoint(filename)

    if level > Fh.RootAttributes['max_level']:
        raise ValueError("Level requested does not exist")

    LDFAB = Fh.SetUpLevelDatasFABs(FABName='data')
    Fh.ReadFABs(LevelDatas=LDFAB, FABName='data')
    compNames=dict([(str(val)[2:-1],int(f[10:])) for f,val in Fh.Root['/'].attrs.items() if 'component_' in f])
    nComp=Fh.Root['/'].attrs['num_components']
    for n in range(nComp):
        compNames[str(Fh.Root['/'].attrs['component_'+str(n)])[2:-1]]=n

    return LDFAB[level],Fh.LevelsAttributes[level]['time'],compNames

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

def main(plotDir,inputFile):

    import os



    try:
        parmDict=getParms(inputname=inputFile)
    except:
        if MPI.COMM_WORLD.Get_rank()==0:
            print(inputFile +" does not exist")
        MPI.COMM_WORLD.Barrier()
        sys.exit(2)

    try:
        os.chdir(plotDir)
    except:
        if MPI.COMM_WORLD.Get_rank()==0:
            print(plotDir +" does not exist")
        MPI.COMM_WORLD.Barrier()
        sys.exit(2)

    from progress.bar import Bar

    CWD = os.getcwd()
    filenames = [f for f in os.listdir(CWD) if os.path.isfile(os.path.join(CWD, f)) and (f[0:4] == 'plot' and '.3d.hdf5' in f)]
    filenames.sort()




    #for key,value in parmDict.items():
    #    print(key," : ", value)
    nu=parmDict['rhs.nu']
    # kappa=parmDict['rhs.bKappa']

    L=parmDict['rhs.L']
    N=parmDict['base.nx']
    Volume=1.
    for l in L:
        Volume*=l
    DX=[l/n for l,n in zip(L,N)]
    perm=((0,1,2),(1,2,0),(2,0,1)) if len(DX)==3 else ((0,1),(1,0))

    KE=np.zeros(len(filenames))
    BF=np.zeros(len(filenames))
    PE=np.zeros(len(filenames))
    diss=np.zeros(len(filenames))
    time=np.zeros(len(filenames))
    with Bar('Processing', max=len(filenames)) as bar:
        for n,filename in enumerate(filenames):


            #try:

            LDFAB,time[n],compNames = GatherComponents(filename)

            LDCoord=LevelDataFAB.MakeFromLDFAB(LDFAB)

            delta=[l/s for l,s in zip(L,N)]
            for x in LDCoord:

                b=x.box
                for p in perm:
                    start=(b.LoEnd()[p[0]])
                    end=start+b.Size()[p[0]]


                    replicate=1
                    for s in p[1:]:
                        replicate*=(b.Size()[s])
                    newshape=tuple([b.Size()[s] for s in p])

                    x[...,p[0]]=np.moveaxis(np.tile((np.arange(start,end)+0.5)*delta[p[0]],replicate).reshape(newshape,order='F'),perm[0],p)


            SpaceDim=LDFAB.SpaceDim
            diss[n]=Diagnostic.epsilon(LDFAB,nu,DX)
            SGSDiss,_,_=Diagnostic.epsilonSGS(LDFAB, LDFAB,DX,nuTComp=compNames['eddyNu'])
            diss[n]+=SGSDiss
            KE[n]=Diagnostic.KEFAB(LDFAB,DX)
            PE[n]=Diagnostic.PE(LDFAB,DX,LDCoord, bComp=compNames['b_total'])
            BF[n]=Diagnostic.BFFAB(LDFAB,DX,wcomp=compNames['z_vel'], bcomp=compNames['b_total'])

            #except:
                #pass
            bar.next()
    for n in range(KE.shape[0]):
        if(MPI.COMM_WORLD.Get_rank()==0):
            if(n>=2):
                budgetKE=abs(((KE[n]-KE[n-2])/(time[n]-time[n-2])+BF[n-1]+diss[n-1])/Volume)
                budgetPE=abs(((PE[n]-PE[n-2])/(time[n]-time[n-2])-BF[n-1])/Volume)
                budgetKE/=max(abs(BF[n-1]),abs(diss[n-1]),abs((KE[n]-KE[n-2])/(time[n]-time[n-2])))/Volume
                budgetPE/=max(abs(BF[n-1]), abs((PE[n]-PE[n-2])/(time[n]-time[n-2])))/Volume
                print("\n" + f'Time =  {time[n-1]: .5e}' +
                                f'; KE =  {KE[n-1]/Volume:.5e}' +
                                f'; PE =  {PE[n-1]/Volume:.5e}' +
                                f'; Diss =  {diss[n-1]/Volume:.5e}'+
                                f'; BF =  {BF[n-1]/Volume:.5e}' +
                                f'; ResidualPE = {budgetPE:.5e}'  +
                                f'; ResidualKE = {budgetKE:.5e}')

    if (MPI.COMM_WORLD.Get_rank()==0):
        integratedBF=integrate.cumtrapz(BF,time)
        integratedDiss=integrate.cumtrapz(diss,time)
        for n in range(0,KE.shape[0]-1):
            print("\n" + f'Time = {time[n]: .5e}' +
                        f'; IntBF+intDiss = {integratedBF[n]+integratedDiss[n]: .5e}' +
                        f'; KE = {KE[n+1]-KE[0]:.5e}')


if __name__=='__main__':
    import getopt
    options="hp:i:"
    long_options=['help', 'plotDir', 'inputFile']
    plotDir=''
    inputFile=''
    try:
        arguments, values = getopt.getopt(sys.argv[1:],options,long_options)
    except getopt.GetoptError:
        if MPI.COMM_WORLD.Get_rank()==0:
            print("Usage: EnergyBudget.py -p <plotDir> -i <inputFile>")
        MPI.COMM_WORLD.Barrier()
        sys.exit(2)

    for argument, value in arguments:
        if argument == '-p':
            plotDir=value
        if argument == '-i':
            inputFile=value
        if argument == '-h':
            if MPI.COMM_WORLD.Get_rank()==0:
                print("Usage: mpirun -np XX EnergyBudget.py -p <plotDir> -i <inputFile> \n")
                print("mpirun -np xx is necessary only if multiple cores are requested")
            MPI.COMM_WORLD.Barrier()
            sys.exit(0)
    main(plotDir, inputFile)
