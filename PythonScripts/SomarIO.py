#!/usr/bin/env python3
import os
from os import path
import sys
import site
sys.path.insert(0,site.USER_SITE)


import h5py as h5
import numpy as np
from FAB import FAB as FAB
from FluxBox import FluxBox
from LevelDataFAB import LevelDataFAB as LDFAB
from LevelDataFB import LevelDataFluxBox as LDFB
from Vects import RealVect, IntVect
from mpi4py import MPI
from Box import Box



class WriteCheckPoint:



    ''' class used to write to and from checkpoints using Chombo data format'''
    def __init__(self, filename):

        try:
            self.FileHandle = h5.File(filename, 'a' if path.exists(filename) else 'w', driver='mpio', comm=MPI.COMM_WORLD, libver='latest')
        except:
            self.FileHandle = h5.File(filename, 'a' if path.exists(filename) else 'w')
        self.Root = self.FileHandle['/']
        self.IsOpen = True


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




        return value

    def WriteLevelDataFAB(self, ldfab, groupname, name, ghosts=(0,0,0), precision='d'):
        # create the dataset
         # first write boxlayout if it does not exist

        q = ldfab

        spaceDim=q.SpaceDim
        ncomp=q.nComp

        if self.IsOpen:

            if groupname not in self.Root.keys(): self.Root.create_group(groupname)
            BLO = [(id, (box.LoEnd() + box.HiEnd())) for id, box in zip(q.PIDs, q.Boxes)]

            if 'boxes' not in self.Root[groupname].keys():
                self.WriteBoxLayout(BLO, groupname)


            if name + "_attributes" not in self.Root[groupname].keys():

                # work on the attributes of the dataset
                self.Root[groupname].create_group(name + "_attributes")
                self.Root[groupname][name + "_attributes"].attrs['comps'] = self.Transmogrify(ncomp)

                hasghosts=ghosts
                outputGhosts=ghosts
                self.Root[groupname][name + "_attributes"].attrs['ghost']=self.Transmogrify(hasghosts)
                self.Root[groupname][name + "_attributes"].attrs['outputGhost']= self.Transmogrify(outputGhosts)
                self.Root[groupname][name + "_attributes"].attrs['objectType'] = self.Transmogrify('FArrayBox')
                #work on creating the datasets for value and offsets


                for key,value in q.SizesForIO.items():
                    self.Root[groupname].create_dataset(name+key,value, dtype=precision)

        # finally each process writes its own fields

        for fab,tag in zip(q,q.tags):
            self.Root[groupname][name+tag][...]=fab[...]


    # writes a fluxbox to the chkpoint file as '/groupname/'
    def WriteLevelDataFluxBox(self,ldvel, groupname, name, ghosts=(1,1,1), precision='d'):
        # use this to start an interactive shell which has access to the local name space
        # from IPython import embed; embed()
        # alternatively use this to start pdb
        # import pdb; pdb.set_trace()
        # note that if running multiple cores, use mpirun -np xx xterm -e "./somarXXX input.file"
        # so that when the embedded shell or the debugger are started, we get a terminal for each instance.

        #try:
        #from IPython import embed; embed()

        if groupname not in self.Root.keys(): self.Root.create_group(groupname)
        vel = ldvel
        SpaceDim = vel.SpaceDim
        ncomp=vel.nComp
        if self.IsOpen:

            if 'boxes' not in self.Root[groupname].keys():
                BLO = [(id, (box.LoEnd() + box.HiEnd())) for id, box in zip(vel.PIDs, vel.Boxes)]

                self.WriteBoxLayout(BLO, groupname)

            # if ncomp != SpaceDim:
            #     raise ValueError("We cannot write FluxBoxes with more than one component per face. This error can be triggered by Output being handled by a process without FABs.")
            if groupname not in self.Root.keys(): self.Root.create_group(groupname)
            if name + "_attributes" not in self.Root[groupname].keys():




                # work on the attributes of the dataset
                self.Root[groupname].create_group(name + "_attributes")
                self.Root[groupname][name + "_attributes"].attrs['comps'] = self.Transmogrify(ncomp)

                hasghosts=ghosts
                outputGhosts=ghosts

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
        #except:
        #    from IPython import embed
        #    embed()











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
        try:
            self.FileHandle = h5.File(filename, 'r', driver='mpio', comm=MPI.COMM_WORLD, libver='latest')
        except:
            self.FileHandle = h5.File(filename, 'r')

        # root
        self.Root = self.FileHandle['/']
        self.RootAttributes = dict(self.FileHandle.attrs.items())
        self.GroupNames = [group for group in self.Root]
        # ChomboGlobal
        self.ChomboGlobalAttributes=dict(self.FileHandle['Chombo_global'].attrs.items())

        self.GenericGroups = [self.Root[group] for group in self.GroupNames if (group[0:6] != 'level_' and group != 'Chombo_global')]
        self.GenericGroupsAttributes = dict()
        # levels
        self.Levels = [self.Root[group] for group in self.GroupNames if group[0:6]=='level_']
        self.LevelsAttributes = list()
        self.LevelsGroupsAttributes = list()
        self.BoxLayout = list()
        #self.DataOffsets = list()
        for group in self.GenericGroups:
            self.GenericGroupsAttributes[group] = (dict(group.attrs.items()))
        for level in self.Levels:
            self.LevelsAttributes.append(dict(level.attrs.items()))
            self.LevelsGroupsAttributes.append(dict([(group, dict(level[group].attrs.items())) for group in level if group[-10:]=='attributes']))
            # boxlayout
            ProcIds = [id % MPI.COMM_WORLD.Get_size() for id in level['Processors'][:]] #
            Boxes = list(level['boxes'][:])
            self.BoxLayout.append([(Id, Box) for (Id,Box) in zip(ProcIds, Boxes)] )
            # get datanames and build dictionary of offsets
            #DataNames = [l[0:-10] for l in level if l[-9:-2] == 'offsets']

            #self.DataOffsets.append(dict([(Name, level[Name + ':offsets=0'][:]) for Name in DataNames]))


    def __str__(self):
        string = "Root attributes \n" + str(self.RootAttributes) + "\n" + \
                 "Chombo_global attributes\n" + str(self.ChomboGlobalAttributes) + "\n" + \
                 "Group Names\n" + str(self.GroupNames) + "\n" + "*****************\n"
        for group in self.GenericGroups:
            string += str(group) + "\n"
            string += str(self.GenericGroupsAttributes[group]) + "\n"
        string += "=======================================================================\n"
        for level in range(len(self.Levels)):
            string += "level " + str(level) + "\n"
            string += "==========================================================================\n"
            string += "LevelsAttributes:\n" + str(self.LevelsAttributes[level]) + "\n\n"
            string += "GroupsAttributes:\n" + str(self.LevelsGroupsAttributes[level]) + "\n\n"
            string += "BoxLayout:\n" + str(self.BoxLayout[level]) + "\n\n"
            string += "Levels keys: \n" + str(self.Levels[level].keys()) + "\n\n"
            #string += "DataOffsets:\n" + str(self.DataOffsets[level])+"\n\n"
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
            ProcId = [(i,) for i, b in self.BoxLayout[n]]


            Domain = tuple(self.LevelsAttributes[n]["domainBox"]) + ('Box',)

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
            ProcId = [(i,) for i, b in self.BoxLayout[n]]

            Domain = tuple(self.LevelsAttributes[n]["domainBox"]) + ('Box',)
            LDFABs = LDFAB((tuple(FABs), tuple(Boxes), tuple(ProcId), Domain,'LevelDataFAB'))
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
            #offsets = self.DataOffsets[n][FABName][:]
            #Data = level[FABName + ":datatype=0"]

            # read the data

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
            # if HasOutputGhosts:
            #     raise ValueError(FABName+ ": we cannot deal yet with outputghosts ")
            ncomps = self.LevelsGroupsAttributes[n][FABName+'_attributes']['comps']
            #offsets = self.DataOffsets[n][FABName][:]
            #Data = level[FABName + ":datatype=0"]

            # read the data
            for fab in LevelDatas[n]:
                b = fab.box

                for (id,valid) in self.BoxLayout[n]:
                    validBox=Box(tuple(valid)+('Box',))
                    if b.ContainsBox(Box(tuple(valid)+('Box',))):
                        if id != comm.Get_rank():
                            raise ValueError("Something is wrong: this box should go to a different process")
                        tag=str(hash(validBox.LoEnd()) )
                        nComps=fab[...].shape[-1]
                        fab[...]=level[FABName+tag][...,0:nComps]



    def IdAndBox(self, i, level):
        return self.BoxLayout[level][i]

    def NumberOfBoxes(self, level):
        return len(self.BoxLayout[level])

if __name__ == "__main__":
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
    Np=(16,8,8)
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





