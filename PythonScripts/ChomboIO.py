
import os
from os import path
import sys
import h5py as h5
import numpy as np
from FAB import FAB as FAB
from LevelDataFAB import LevelDataFAB as LDFAB
from LevelDataFB import LevelDataFluxBox as LDFB
from mpi4py import MPI
from Box import Box
import PyGlue as pg
from SOMAR import SOMAR

from Vects import IntVect, RealVect
class WriteCheckPoint:
    ''' class used to write to and from checkpoints (but also plot files) using Chombo data format'''
    def __init__(self, filename):
        # do NOT use libver='latest'. Causes issues with Paraview
        if True:
            self.FileHandle = h5.File(filename, 'w', driver='mpio', comm=MPI.COMM_WORLD)
            self.Root = self.FileHandle['/']
            self.IsOpen = True
        else:
            self.IsOpen = False


    def __del__(self):
        if self.IsOpen:
            self.FileHandle.close()

    def WriteHeader(self, Group, key, value):
        if not self.IsOpen: return

        if Group[0] == '/':

            self.Root.attrs[key] = self.Transmogrify(value)
        else:

            if Group not in self.Root.keys(): self.Root.create_group(Group)
            self.Root[Group].attrs[key] = self.Transmogrify(value)


    def Transmogrify(self, value):
        if isinstance(value, str):
            return np.bytes_(value)

        if isinstance(value, int):
            return np.int32(value)

        if isinstance(value, Box):
            Dim=len(value.LoEnd())
            if Dim == 2:
                names = ['lo_i', 'lo_j', 'hi_i', 'hi_j']
            else:
                names = ['lo_i', 'lo_j', 'lo_k', 'hi_i', 'hi_j', 'hi_k']
            dt = np.dtype({'names': names,
                        'formats': 2*Dim*['<i4'],
                        'offsets': [4*i for i in range(2*Dim)],
                        'itemsize': (2*Dim + 1)* 4})

            return np.array(value.LoEnd()+value.HiEnd(), dtype=dt)
        
        if isinstance(value, IntVect):
            names = ['intvecti', 'intvectj']
            dim=len(value.Vect)
            if dim ==3: names.append('intvectk')
            
            dt = np.dtype({'names': names,
                        'formats': dim*['<i4'],
                        'offsets': [4*i for i in range(dim)],
                        'itemsize': (dim + 1)* 4})
            return np.array(value.Vect, dtype=dt)
        
        if isinstance(value, RealVect):
            names = ['x', 'y']
            dim=len(value.Vect)
            if dim ==3: names.append('z')
            
            dt = np.dtype({'names': names,
                        'formats': dim*['<f8'],
                        'offsets': [8*i for i in range(dim)],
                        'itemsize': (dim + 1)* 8})
            return np.array(value.Vect, dtype=dt)
        # try:

            


        #     if value[-1] == "IntVect":
        #         names = ['intvecti', 'intvectj']
        #         if len(value[:-1]) == 3: names.append('intvectk')
        #         dt = np.dtype({'names': names,
        #                 'formats': len(value[:-1])*['<i4'],
        #                 'offsets': [4*i for i in range(len(value[:-1]))],
        #                 'itemsize': (len(value[:-1]) + 1)* 4})
        #         return np.array(value[:-1], dtype=dt)

        #     if value[-1] == "RealVect":
        #         names = ['x', 'y']
        #         if len(value[:-1]) == 3: names.append('z')
        #         dt = np.dtype({'names': names,
        #                 'formats': len(value[:-1])*['<f8'],
        #                 'offsets': [8*i for i in range(len(value[:-1]))],
        #                 'itemsize': (len(value[:-1]) + 1)* 8})

        #         return np.array(value[:-1], dtype=dt)
        # except:
        return value

    def WriteLevelDataFAB(self, ldfab, groupname, name, ghosts=(0,0,0), precision='d'):
        # create the dataset
        # first write boxlayout if it does not exist
        # Note that we rely on root to broadcast commom parameters. 
        # this to avoid potential issues with IO::WRITEHDF5
        
        q = ldfab
        
        
        if self.IsOpen:
            # if ncomp is 0:
            #     raise ValueError(" root has not fabs. Output should be assigned to a process with fabs")
            if groupname not in self.Root.keys(): self.Root.create_group(groupname)
            
            if 'boxes' not in self.Root[groupname].keys():
                if MPI.COMM_WORLD.Get_rank()==0:
                    BLO = [(id, (box.LoEnd() + box.HiEnd())) for id, box in zip(q.PIDs, q.Boxes)]
                else:
                    BLO=None
                BLO=MPI.COMM_WORLD.bcast(BLO, root=0)

                self.WriteBoxLayout(BLO, groupname)

            if MPI.COMM_WORLD.Get_rank()==0:
                    size = 0
                    offsets=list()
                    for b in q.Boxes:
                        box = Box.FromAnotherBox(b)
                        box.Grow(ghosts)
                        offsets.append(size)
                        size += q.nComp * box.Ncells()

                    offsets.append(size)
            else:
                offsets=None

            offsets=MPI.COMM_WORLD.bcast(offsets,root=0)
            size=offsets[-1]
            
            if name not in self.Root[groupname].keys():
                if MPI.COMM_WORLD.Get_rank()==0:
                    b = [box for box in q.Boxes if q.FABs[0].box.ContainsBox(box)][0]
                    hasghosts= IntVect(tuple([int(s2>s1) for s1, s2 in zip(b.Size(), q.FABs[0].box.Size())])+('IntVect',))
                    outputGhosts = IntVect(tuple([ghosts[i] for i in range(len(b.Size()))])+('IntVect',))
                    if any([ g<oG for g, oG in zip(ghosts, outputGhosts)]):
                        print(ghosts, outputGhosts)
                        raise ValueError(" FAB without ghosts is being written with outputGhosts")
                else:
                    hasghosts,outputGhosts=None,None

                hasghosts=MPI.COMM_WORLD.bcast(hasghosts,root=0)
                outputGhosts=MPI.COMM_WORLD.bcast(outputGhosts,root=0)
                
                # work on the attributes of the dataset
                self.Root[groupname].create_group(name + "_attributes")
                self.Root[groupname][name + "_attributes"].attrs['comps'] = self.Transmogrify(q.nComp)
                self.Root[groupname][name + "_attributes"].attrs['ghost']=self.Transmogrify(hasghosts)
                self.Root[groupname][name + "_attributes"].attrs['outputGhost']= self.Transmogrify(outputGhosts)
                self.Root[groupname][name + "_attributes"].attrs['objectType'] = self.Transmogrify('FArrayBox')

                #work on creating the datasets for value and offsets
                self.Root[groupname].create_dataset(name + ":datatype=0", (size,), dtype=precision)
                self.Root[groupname].create_dataset(name + ":offsets=0", data=offsets)

        
            for i, (b,offset) in enumerate(zip(q.Boxes,offsets)):

                box = Box.FromAnotherBox(b)
                box.Grow(ghosts)
                 
                for fab in q:
                    if fab.box.ContainsBox(box):
                        try:
                            
                            size=q.nComp*fab.box.Ncells()
                            self.Root[groupname][name + ":datatype=0"][offset:offset+size] = fab.data.ravel(order='F')
                        except:
                            #pass
                            dummy=FAB.FromBox(box,ncomp=q.nComp)
                            dummy.copyFrom(fab,destComp=tuple(range(q.nComp)), srcComp=tuple(range(q.nComp)))
                            size=q.nComp*box.Ncells()
                            self.Root[groupname][name + ":datatype=0"][offset:offset+size] = dummy.data.ravel(order='F')

             




    # writes a fluxbox to the chkpoint file as '/groupname/'
    def WriteLevelDataFluxBox(self,ldvel, groupname, name, ghosts=(0,0,0), precision='d'):

        vel = ldvel
        SpaceDim = vel.SpaceDim
        if len(vel.FBs) >0:
            ncomp = vel.FBs[0].Fluxes[0].data.shape[-1::][0] * SpaceDim
        else:
            ncomp = 0  # this should trigger an error below if it ever happens that root has not boxes
        if self.IsOpen: #only root writes.

            if 'boxes' not in self.Root[groupname].keys():
                BLO = [(id, (box.LoEnd() + box.HiEnd())) for id, box in zip(vel.PIDs, vel.Boxes)]
                self.WriteBoxLayout(BLO, groupname)

            if ncomp != SpaceDim:
                raise ValueError("We cannot write FluxBoxes with more than one component per face. This error can be triggered by Output being handled by a process without FABs.")
            if groupname not in self.Root.keys(): self.Root.create_group(groupname)
            if name not in self.Root[groupname].keys():

                size = 0
                offsets=list()
                for b in vel.Boxes:
                    offsets.append(size)
                    size += ncomp * b.Ncells()

                offsets.append(size)


                # work on the attributes of the dataset
                self.Root[groupname].create_group(name + "_attributes")
                self.Root[groupname][name + "_attributes"].attrs['comps'] = self.Transmogrify(ncomp)
                b = [box for box in vel.Boxes if vel.FBs[0].box.ContainsBox(box)][0]
                hasghosts=IntVect(tuple([s2>s1 for s1, s2 in zip(b.Size(), vel.FBs[0].box.Size())])+('IntVect',))
                outputGhosts = IntVect(tuple([ghosts[i] for i in range(len(b.Size()))])+('IntVect',))
                self.Root[groupname][name + "_attributes"].attrs['ghost'] = self.Transmogrify(hasghosts)
                self.Root[groupname][name + "_attributes"].attrs['outputGhost'] = self.Transmogrify(outputGhosts)
                self.Root[groupname][name + "_attributes"].attrs['objectType'] = self.Transmogrify('unknown')

                #work on creating the datasets for value and offsets


                for key,value in vel.SizesForIO.items():
                    self.Root[groupname].create_dataset(name+key,value, dtype=precision)




                for v,tags in zip(vel, vel.tags):
                    for flux,tag in zip(v.Fluxes,tags):
                        self.Root[groupname][name+tag][...]=flux[...]


        MPI.COMM_WORLD.barrier()
        #         self.Root[groupname].create_dataset(name + ":datatype=0", (size,), dtype=precision)
        #         self.Root[groupname].create_dataset(name + ":offsets=0", data=offsets)


        # for i,b in enumerate(vel.Boxes):
        #     X=vel.GatherFluxBox(b)
        #     if self.IsOpen:
        #         begin=self.Root[groupname][name+":offsets=0"][i]
        #         for fab in X.Fluxes:
        #             S=1
        #             for s in fab.data.shape:
        #                 S*=s
        #             end=begin+S
        #             self.Root[groupname][name+":datatype=0"][begin:end] = fab.data.ravel(order='F')
        #             begin = end




    def WriteBoxLayout(self, BLO, group):

        if self.IsOpen:
            SpaceDim = len(BLO[0][1])//2

            if SpaceDim == 2:
                names = ['lo_i', 'lo_j', 'hi_i', 'hi_j']
            else:
                names = ['lo_i', 'lo_j', 'lo_k', 'hi_i', 'hi_j', 'hi_k']

            dt = np.dtype({'names': names,
                           'formats': (2*SpaceDim)*['<i4'],
                           'offsets': [4*i for i in range(2*SpaceDim)],
                           'itemsize': (2 * SpaceDim + 1) * 4})

            data=np.array([tuple(box) for (id, box) in BLO], dtype=dt)
            self.Root[group].create_dataset('boxes', data=data)
            self.Root[group].create_dataset('Processors', data=[id for (id, box) in BLO])

class ReadCheckPoint:
    """ handles reading checkpoint files, providing SOMAR with info to buid LevelData<T> objects stored in the file"""
    def __init__(self, filename):
        self.FileHandle = h5.File(filename, 'r', driver='mpio', comm=MPI.COMM_WORLD)
        # root
        self.Root = self.FileHandle['/']
        self.RootAttributes = dict([(attribute, self.FileHandle.attrs[attribute]) for attribute in self.FileHandle.attrs])
        self.GroupNames = [group for group in self.Root]
        # ChomboGlobal
        self.ChomboGlobalAttributes=dict([(attribute, self.FileHandle[self.GroupNames[0]].attrs[attribute]) for attribute in self.FileHandle['Chombo_global'].attrs])
        # levels

        self.Levels = [self.Root[group] for group in self.GroupNames if group[0:6]=='level_']
        self.LevelsAttributes = list()
        self.LevelsGroupsAttributes = list()
        self.BoxLayout = list()
        self.DataOffsets = list()
        for level in self.Levels:
            self.LevelsAttributes.append(dict([(attribute, level.attrs[attribute]) for attribute in level.attrs]))
            self.LevelsGroupsAttributes.append(dict([(group, dict([(attribute, level[group].attrs[attribute]) for attribute in level[group].attrs])) for group in level if group[-10:]=='attributes']))
            # boxlayout
            ProcIds = [id % MPI.COMM_WORLD.Get_size() for id in level['Processors'][:]] #
            Boxes = list(level['boxes'][:])
            self.BoxLayout.append([(Id, Box) for (Id,Box) in zip(ProcIds, Boxes)] )
            # get datanames and build dictionary of offsets
            DataNames = [l[0:-10] for l in level if l[-9:-2] == 'offsets']

            self.DataOffsets.append(dict([(Name, level[Name + ':offsets=0'][:]) for Name in DataNames]))

    def __str__(self):
        string = "Root attributes \n" + str(self.RootAttributes) + "\n" + \
                 "Chombo_global attributes\n" + str(self.ChomboGlobalAttributes) + "\n" + \
                 "Group Names\n" + str(self.GroupNames) + "\n" + "*****************\n"

        for level in range(len(self.Levels)):
            string += "level " + str(level) + "\n"
            string += "==========================================================================\n"
            string += "LevelsAttributes:\n" + str(self.LevelsAttributes[level]) + "\n\n"
            string += "GroupsAttributes:\n" + str(self.LevelsGroupsAttributes[level]) + "\n\n"
            string += "BoxLayout:\n" + str(self.BoxLayout[level]) + "\n\n"
            string += "DataOffsets:\n" + str(self.DataOffsets[level])+"\n\n"
            string += "==========================================================================\n"

        return string

    def __del__(self):
        self.FileHandle.close()


    def SetUpLevelDatasFluxBoxes(self, FABName='velData'):

        comm = MPI.COMM_WORLD
        LevelDatas = list()
        for n in range(len(self.Levels)):
            FBs = list()
            IsFB = str(self.LevelsGroupsAttributes[n][FABName + '_attributes']['objectType'])[2:-1] == 'unknown'
            if not IsFB:
                raise ValueError(FABName+" is not a fluxbox")
            HasGhosts = self.LevelsGroupsAttributes[n][FABName+'_attributes']['ghost']
            ncomps = self.LevelsGroupsAttributes[n][FABName+'_attributes']['comps']


            for (id, box) in self.BoxLayout[n]:
                if id == comm.Get_rank():
                    # determine if box needs ghosts
                    Args = []


                    b = Box(tuple(box)+('Box',))
                    LoEnd = []
                    Size = []
                    for lo, size, g in zip(b.LoEnd(), b.Size(), HasGhosts):
                        LoEnd.append(lo - g)
                        Size.append(size + 2*g)
                    # pack arguments
                    Args.append(tuple(LoEnd))
                    Args.append(tuple(Size))
                    Args.append(b.Centering())
                    Args.append((ncomps,))
                    Args.append('FluxBox')
                    FBs.append(Args)

            Boxes = [tuple(b)+('Box',) for i, b in self.BoxLayout[n]]
            ProcId=[(i,) for i,b in self.BoxLayout[n]]
            Domain=tuple(self.LevelsAttributes[n]['domainBox'])+('Box',)
            LevelDatas.append(LDFB((tuple(FBs), tuple(Boxes), tuple(ProcId), Domain,'LevelDataFluxBox')))

        return LevelDatas

    def SetUpLevelDatasFABs(self, FABName='qData'):
        comm = MPI.COMM_WORLD
        LevelDatas = list()
        for n in range(len(self.Levels)):
                # get # of processors



            # define the local FABs needed
            FABs = list()
            IsFAB = str(self.LevelsGroupsAttributes[n][FABName + '_attributes']['objectType'])[2:-1] == 'FArrayBox'
            if not IsFAB:
                raise ValueError(FABName+" is not a FArrayBox")
            HasGhosts = self.LevelsGroupsAttributes[n][FABName+'_attributes']['ghost']
            ncomps = self.LevelsGroupsAttributes[n][FABName+'_attributes']['comps']


            for (id, box) in self.BoxLayout[n]:
                if id == comm.Get_rank():
                    # determine if box needs ghosts
                    Args = []
                    b = Box(tuple(box)+('Box',))
                    LoEnd = []
                    Size = []
                    for lo, size, g in zip(b.LoEnd(), b.Size(), HasGhosts):
                        LoEnd.append(lo - g)
                        Size.append(size + 2*g)
                    # pack arguments
                    Args.append(tuple(LoEnd))
                    Args.append(tuple(Size))
                    Args.append(b.Centering())
                    Args.append((ncomps,))
                    Args.append(None)
                    Args.append('FAB')
                    FABs.append(Args)

            Boxes = [tuple(b)+('Box',) for i, b in self.BoxLayout[n]]
            ProcId=[(i,) for i,b in self.BoxLayout[n]]

            Domain=tuple(self.LevelsAttributes[n]['prob_domain'])+('Box',)

            LDFABs = LDFAB((tuple(FABs), tuple(Boxes), tuple(ProcId), Domain, 'LevelDataFAB'))
            LevelDatas.append(LDFABs)

        return LevelDatas

    def ReadFluxBoxes(self, LevelDatas=None, FABName='velData'):
        comm = MPI.COMM_WORLD

        SpaceDim = self.ChomboGlobalAttributes['SpaceDim']
        if LevelDatas==None:
            LevelDatas = self.SetUpLevelDatasFluxBoxes(FABName)
            self.ReadFABs(LevelDatas=LevelDatas, FABName=FABName)

        for n,level in enumerate(self.Levels):

            IsFB = str(self.LevelsGroupsAttributes[n][FABName + '_attributes']['objectType'])[2:-1] == 'unknown'
            if not IsFB:
                raise ValueError(FABName+" is not a Fluxbox")
            HasOutputGhosts = self.LevelsGroupsAttributes[n][FABName+'_attributes']['outputGhost']
            ncomps = self.LevelsGroupsAttributes[n][FABName+'_attributes']['comps']
            # offsets = self.DataOffsets[n][FABName][:]
            # Data = level[FABName + ":datatype=0"]

            # read the data
            # for fluxbox in LevelDatas[n]:

            #     b = fluxbox.box

            #     for i, (id, valid) in enumerate(self.BoxLayout[n]):

            #         ofst = offsets[i]
            #         if b.ContainsBox(Box(tuple(valid)+('Box',))):
            #             if id != comm.Get_rank():
            #                 raise ValueError("Something is wrong: this box should go to a different process")
            #             for d in range(SpaceDim):
            #                 fab = fluxbox.Fluxes[d]
            #                 if SpaceDim == 3:
            #                     size = (fab.data.shape[0] -2) * (fab.data.shape[1] - 2) * (fab.data.shape[2] - 2)
            #                     S = (fab.data.shape[0]-2, fab.data.shape[1] - 2, fab.data.shape[2] - 2)
            #                     fab[1:-1, 1:-1, 1:-1, 0] = Data[ofst:ofst + size].reshape(S, order='F')
            #                     ofst+=size
            #                 else:
            #                     size = (fab.data.shape[0] -2) * (fab.data.shape[1] - 2)
            #                     S = (fab.data.shape[0]-2, fab.data.shape[1] - 2)
            #                     fab[1:-1, 1:-1, 0] = Data[ofst:ofst + size].reshape(S, order='F')
            #                     ofst += size
            #             if ofst > offsets[i + 1]:
            #                 print(ofst, " ", offsets[i+1], " ", (offsets[i+1]-offsets[i]), " ", size , ' ', fab.data.shape )
            #                 raise ValueError("mismatch between data used and offsets.")


            for fluxbox in LevelDatas[n]:

                b = fluxbox.box

                for  (id, valid) in self.BoxLayout[n]:
                    validBox=Box(tuple(valid)+('Box',))
                    #ofst = offsets[i]
                    if b.ContainsBox(Box(tuple(valid)+('Box',))):
                        if id != comm.Get_rank():
                            raise ValueError("Something is wrong: this box should go to a different process")
                        for d in range(SpaceDim):
                            fab = fluxbox.Fluxes[d]
                            tag=str(hash(validBox.LoEnd()+(d,)) )

                            fab[...]=level[FABName+tag][...]

    def ReadFABs(self,LevelDatas=None,FABName='qData'):

        comm = MPI.COMM_WORLD

        SpaceDim = self.ChomboGlobalAttributes['SpaceDim']
        if LevelDatas==None:
            LevelDatas = self.SetUpLevelDatasFABs(FABName)
            self.ReadFABs(LevelDatas=LevelDatas, FABName=FABName)

        for n,level in enumerate(self.Levels):

            IsFAB = str(self.LevelsGroupsAttributes[n][FABName + '_attributes']['objectType'])[2:-1] == 'FArrayBox'
            if not IsFAB:
                raise ValueError(FABName+" is not a FArrayBox")
            HasOutputGhosts = any(tuple(self.LevelsGroupsAttributes[n][FABName+'_attributes']['outputGhost']))
            #if HasOutputGhosts:
            #    raise ValueError(FABName+ ": we cannot deal yet with outputghosts ")
            ncomps = self.LevelsGroupsAttributes[n][FABName+'_attributes']['comps']
            offsets = self.DataOffsets[n][FABName][:]
            Data = level[FABName + ":datatype=0"]

            # read the data
            for fab in LevelDatas[n]:
                b = fab.box

                for i, (id,valid) in enumerate(self.BoxLayout[n]):
                    if b.ContainsBox(Box(tuple(valid)+('Box',))):
                        if id != comm.Get_rank():
                            raise ValueError("Something is wrong: this box should go to a different process")
                        if HasOutputGhosts:
                            fab[...] = Data[offsets[i]:offsets[i + 1]].reshape(tuple(b.Size())+(ncomps,), order='F')
                        else:
                            if SpaceDim==3:
                                fab[1:-1, 1:-1, 1:-1,:] = Data[offsets[i]:offsets[i + 1]].reshape(tuple(Box(tuple(valid)+('Box',)).Size()) + (ncomps,), order='F')
                            else:
                                fab[1:-1, 1:-1,:] = Data[offsets[i]:offsets[i + 1]].reshape(tuple(Box(tuple(valid)+('Box',)).Size()) + (ncomps,), order='F')

    def IdAndBox(self, i, level):
        return self.BoxLayout[level][i]

    def NumberOfBoxes(self, level):
        return len(self.BoxLayout[level])


def main():
    if MPI.COMM_WORLD.Get_rank()==0:
        if path.exists('test.hdf5'):
            os.remove('test.hdf5')
    MPI.COMM_WORLD.barrier()
    # create a leveldata
    # create a set of boxes and Proc Ids, the latter going round robin
    valid=[]
    Id=[]
    FB=[]
    # the size of the domain is (Np[0]*boxSize[0], Np[1]*boxSize[1],...)
    Np=(1,1,1)
    boxSize=(16,16,32)
    for k in range(Np[2]):
        for j in range(Np[1]):
            for i in range(Np[0]):
                loEnd=(i*boxSize[0], j*boxSize[1], k*boxSize[2])
                hiEnd=((i+1)*boxSize[0]-1, (j+1)*boxSize[1]-1,(k+1)*boxSize[2]-1)
                valid.append(loEnd+hiEnd+("Box",))
                Id.append(((i+j*Np[0]+k*Np[0]*Np[1]) % MPI.COMM_WORLD.Get_size(),))

                if Id[-1][0]==MPI.COMM_WORLD.Get_rank(): # create a Flux Box
                    Args=[]
                    Args.append((loEnd[0]-1,loEnd[1]-1,loEnd[2]-1))
                    Args.append((boxSize[0]+2,boxSize[1]+2,boxSize[2]+2))
                    Args.append((0,0,0))
                    Args.append((1,))
                    Args.append('FluxBox')
                    FB.append(Args)
    Domain = (0,0,0)+(Np[0]*boxSize[0]-1,Np[1]*boxSize[1]-1,Np[2]*boxSize[2]-1) + ('Box',)
    ldFluxBoxOut=LDFB((FB,valid,Id,Domain,"LevelDataFluxBox"))

    for fb in ldFluxBoxOut:
        for d in range(3):
            fb.Fluxes[d][...]=float(MPI.COMM_WORLD.Get_rank()+d)


    W=WriteCheckPoint('test.hdf5')

    W.WriteHeader('Chombo_global','SpaceDim',3)
    W.WriteHeader('Chombo_global', 'testReal',float(0))
    W.WriteHeader('level_0','dt',0.)
    W.WriteHeader('level_0','domainBox', Box(Domain))
    W.WriteLevelDataFluxBox(ldFluxBoxOut,'level_0','velData')
    del W
    MPI.COMM_WORLD.barrier()
    R=ReadCheckPoint('test.hdf5')
    if MPI.COMM_WORLD.Get_rank()==0: print(R)
    #if MPI.COMM_WORLD.Get_rank()==0: print(R)
    ldFluxBoxIn=R.SetUpLevelDatasFluxBoxes()
    R.ReadFluxBoxes(LevelDatas=ldFluxBoxIn)
    errors=[]
    for (vout,vin) in zip(ldFluxBoxOut,ldFluxBoxIn[0]):

        for d in range(3):
            errors.append(np.max(abs(vin.Fluxes[d][...]-vout.Fluxes[d][...])))

    if MPI.COMM_WORLD.Get_rank()>0:
        MPI.COMM_WORLD.send(errors,dest=0, tag=1)
    else:
        print("Proc 0", "errors ", errors)
        for k in range(MPI.COMM_WORLD.Get_size()-1):
            dummy=MPI.COMM_WORLD.recv(source=k+1, tag=1)
            print("Proc ", k+1, "errors ", dummy)

    MPI.COMM_WORLD.Barrier()
    del R
    MPI.COMM_WORLD.Barrier()
    if MPI.COMM_WORLD.Get_rank()==0:
        if path.exists('test.hdf5'):
            os.remove('test.hdf5')

if __name__ == "__main__":
    main()

