
import os
from os import path
from LevelDataFAB import LevelDataFAB
from Box import Box



from SOMAR import SOMAR
import PyPhysics
import Slicer
# check for paths and if not make them

NPYDIR="npy"

NPYDIR=NPYDIR+"/"

#if not os.path.exists(NPYDIR):
#    os.mkdir(NPYDIR)

counter = 0
parmDict = PyPhysics.Physics().GetParameters()
comm = PyPhysics.Physics().GetComm()
MpiNumProcs = PyPhysics.Physics().GetNumProcs()
myMpiId = comm.Get_rank()


@SOMAR
def PostStep(vel, p, q, time, dir, pos):
    global counter
    counter+=1
    output_freq=1 # modify this value to suit
    if counter%output_freq==0:

        Q=Slicer.GatherSlice(q, dir=dir, pos=pos)
        P=Slicer.GatherSlice(p, dir=dir, pos=pos)
        V=Slicer.GatherSlice(vel, dir=dir, pos=pos)

        if myMpiId==0:
            filename=NPYDIR+'Q_iter_'+str(counter).zfill(7)+'_dir_'+str(dir)+'_pos_'+str(pos)+'.npy'
            try:
                Q.Write(filename)
            except: pass
            finally:
                del Q

            filename=NPYDIR+'P_iter_'+str(counter).zfill(7)+'_dir_'+str(dir)+'_pos_'+str(pos)+'.npy'
            try:
                P.Write(filename)
            except: pass
            finally:
                del P

            filename=NPYDIR+'V_iter_'+str(counter).zfill(7)+'_dir_'+str(dir)+'_pos_'+str(pos)+'.npz'
            try:
                V.Write(filename)
            except: pass
            finally:
                del V




