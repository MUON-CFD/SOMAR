

import os
from os import path
import sys
import numpy as np
from FAB import FAB as FAB
from LevelDataFAB import LevelDataFAB as LDFAB
from LevelDataFB import LevelDataFluxBox as LDFB
from Vects import IntVect 
from mpi4py import MPI
from Box import Box
import PyGlue as pg
from SOMAR import SOMAR
from BoxLayout import BoxLayout as BL
try:
    import matplotlib.pyplot as plt
except:
    pass

import SomarIO # import this one first!
import ChomboIO
WCP = SomarIO.WriteCheckPoint

RCP = SomarIO.ReadCheckPoint

WPF = ChomboIO.WriteCheckPoint

FileHandle = dict()

@SOMAR
def OpenFileForWrite(filename, checkpoint=False):
    
    
    if checkpoint:

        FileHandle[filename]=WCP(filename)
    else:

        FileHandle[filename]=WPF(filename)
    
    

@SOMAR
def CloseFile(filename):
    del FileHandle[filename]

@SOMAR
def WriteToHeader(filename,groupname, key, value):
    if filename in FileHandle:
        FileHandle[filename].WriteHeader(groupname, key, value)
    else:
        raise ValueError("attempting to write to Header of file " + filename + "but file needs to opened first")
        

# @SOMAR
# def OpenFileForRead(filename):

#     FileHandle[filename]=RCP(filename)

@SOMAR
def ReadFromHeader(filename,groupname, Key):
    
    if groupname[0]=="/":
        S= FileHandle[filename].RootAttributes[Key]
    else:

        S = FileHandle[filename].Root[groupname].attrs[Key]
    # need to do the following otherwise tuples are not upacked properly.
    R = S
    
     # need to do the following otherwise tuples are not upacked properly.
    
    if (type(S) is np.ndarray ):
        R=tuple(S)
    if (type(S) is np.void):
        R=tuple(S)
    if (type(R) is np.bytes_):
        R = str(R)[2:-1]
    
    
    return R



def FileDump(filename):
    print(FileHandle[filename])

def GetHeaderString(filename):
    return FileHandle[filename].__str__()

@SOMAR
def OpenFileForRead(filename):
    try:

        FileHandle[filename]=RCP(filename)
        return 0
    except:
        return -1
@SOMAR
def SizeOfBoxLayout(filename,level):
    R = FileHandle[filename].BoxLayout[level]
    
    return len(R)

@SOMAR
def GetBox(filename,level, position):
    
    R = FileHandle[filename].BoxLayout[level]
    
    return tuple([x.item() for x in R[position][1]])
@SOMAR
def getProcID(filename,level, position):
    R = FileHandle[filename].BoxLayout[level]
    return R[position][0]

@SOMAR
def ReadFAB(filename,level, LD, name):
    R = FileHandle[filename]
    LevelsDataFAB=R.SetUpLevelDatasFABs(FABName=name)
    R.ReadFABs(LevelDatas=LevelsDataFAB, FABName=name)


    for localfab in LevelsDataFAB[level]:
        for fab in LD:
            if localfab.box.ContainsBox(fab.box):

                fab.copyFrom(localfab, destComp=list(range(fab.ncomp)), srcComp=list(range(localfab.ncomp)))


@SOMAR
def ReadFB(filename,level, LD, name):
    R = FileHandle[filename]

    LevelsDataFB = R.SetUpLevelDatasFluxBoxes(FABName=name)
    R.ReadFluxBoxes(LevelDatas=LevelsDataFB, FABName=name)

    for localfb in LevelsDataFB[level]:
        for fb in LD:
            if localfb.box.ContainsBox(fb.box):
                fb.copyFrom(localfb)


@SOMAR
def WriteLevelDataFAB(filename, groupname,  name, ldfab, ghosts):
    
    
    FileHandle[filename].WriteLevelDataFAB( ldfab, groupname,name,ghosts=ghosts)
    
    

@SOMAR
def WriteCheckPointLevelDataFAB(filename,groupname,name,ldfab,ghosts):
    
    FileHandle[filename].WriteLevelDataFAB( ldfab, groupname,name,ghosts=ghosts)
    

@SOMAR
def WriteCheckPointLevelDataFluxBox(filename, groupname, name, ldfb, ghosts):

    
    FileHandle[filename].WriteLevelDataFluxBox(ldfb, groupname, name, ghosts=ghosts)
    

@SOMAR
def DeleteIfFileExists(filename):
    if MPI.COMM_WORLD.Get_rank()==0:
        if path.exists(filename):
            os.remove(filename)
    MPI.COMM_WORLD.barrier()

@SOMAR
def WriteBLToFile(filename,bl):
    
    bl.write(filename)

@SOMAR 
def ReadBLFromFile(filename):
    return BL.read(filename).toCPP()

@SOMAR
def printBL(bl):
    
    print(bl)
    bl.write('test')
    
def PrintBL(bl):

    print(bl)
    bl.write('test')

def main(howmany=1):
    
    from random import seed
    from random import randint
    if MPI.COMM_WORLD.Get_rank()==0:
        if path.exists('test.hdf5'):
            os.remove('test.hdf5')
    # seed random number generator
    seed(1)
    count=0
    nx, ny, nz = (256, 128, 512) # total domain
    mx, my, mz = (256//4, 128//8, 512//32)  # min block
    while(count<howmany):
        # create collection of boxes
        comm = MPI.COMM_WORLD
        MyRank = comm.Get_rank()
        MySize = comm.Get_size()

        Boxes = []
        Pids = []
        FABs = []
        # Generate the boxes that cover the domain and assign them to a random process
        for k in range(nz // mz):
            for j in range(ny // my):
                for i in range(nx // mx):
                    b = Box.FromLoEndAndSize((i * mx, j * my, k * mz), (mx, my, mz))
                    Boxes.append(b)

                    Pids.append((randint(0, MySize - 1), ))


        # create a fabs
        for b, id in zip(Boxes, Pids):

            if id[0] == MyRank:
                Args = []

                G = []
                S = []
                for l, s in zip(b.LoEnd(), b.Size()):
                    G.append(l - 1)
                    S.append(s + 2)

                Args.append(tuple(G))
                Args.append(tuple(S))
                Args.append((0, 0, 0))
                Args.append((1,))
                Args.append(None)
                Args.append('FAB')
                FABs.append(Args)

        Args=5*[None]
        Args[4] = 'LevelDataFAB'
        Args[3] = (0,0,0, nx-1,ny-1,nz-1, 'Box')
        Args[2] = Pids
        Args[1] = [(b.LoEnd() + b.HiEnd()+('Box',)) for b in Boxes]
        Args[0] = FABs

        LD = LDFAB(Args)

        OpenFileForWrite('test.hdf5')
        WriteToHeader("level_0",'test',0)
        CloseFile()
        WriteLevelDataFAB('test.hdf5', "level_0", "pData", LD, IntVect((1,1,1)+('IntVect',)))
        if MPI.COMM_WORLD.Get_rank()==0:
            if path.exists('test.hdf5'):
                os.remove('test.hdf5')
        count+=1



    pass

if __name__ == "__main__":
    main(100)
