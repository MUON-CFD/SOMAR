import numpy as np
import PyGlue as pg
from Box import Box
from FAB import FAB
from mpi4py import MPI
import h5py as h5
from LevelDataFAB import LevelDataFAB as LD
from LevelDataFB import LevelDataFluxBox as LDFB
import matplotlib.pyplot as plt
# test reading  checkpoint file
# hierarchy of file
# GROUP "/" {
#  GROUP "Chombo_global",
#  GROUP "level_0"{
#         GROUP "pData_attributes",
#         GROUP "qData_attributes",
#         GROUP "velData_attributes"
#        }
#    GROUP "level_1"{
#           GROUP "pData_attributes",
#           GROUP "qData_attributes",
#           GROUP "velData_attributes"
#         }
#    etc.
# }
# RootAttributes = ['base.blockFactor', 'time.dtMult',
#                     'time.fixedDt', 'base.isPeriodic_0',
#                     'base.isPeriodic_1', 'base.isPeriodic_2',
#                     'time.maxDt', 'amr.maxLevel',
#                     'base.nx', 'base.nxOffset',
#                     'amr.useSubcycling',
#                     'iteration','max_level',
#                     'numQComps', 'numScalars',
#                     'num_levels', 'regrid_interval_0',
#                     'rhs.L', 'rhs.bKappa',
#                     'rhs.coriolisF', 'rhs.eddyPrandtl',
#                     'rhs.eddyViscMethodSize', 'rhs.eddyViscMethod_0',
#                     'rhs.nu', 'rhs.sKappaSize',
#                     'rhs.sKappa_0', 'scalarComponent_0',
#                     'time']
#ChomboGlobalAttributes = ['SpaceDim', 'testReal']
# LevelsAttributes = ['domainBox', 'dt',
#                         'finestExtantLevel', 'isEmpty',
#                         'ref_ratio', 'time', 'vec_dx']

# LevelDatasets = ['Processors', 'boxes',
#                 'pData:datatype=0', 'pData:offsets=0',
#                 'qData:datatype=0',
#                 'qData:offsets=0',
#                 'velData:datatype=0', 'velData:offsets=0',
#                 ]

# LevelGroups = ['pData_attributes', 'qData_attributes', 'velData_attributes']

# LevelGroupsAttributes = ['comps', 'ghost', 'objectType', 'outputGhost']
class WriteCheckPoint:
    def __init__(self, filename):

        if MPI.COMM_WORLD.Get_rank()==0:
            self.FileHandle = h5.File(filename, 'a')
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
            self.Root.attrs[key] = value
        else:
            if Group not in self.Root.keys(): self.Root.create_group(Group)
            self.Root[Group].attrs[key] = value




    def WriteLevelDataFAB(self, ldfab, groupname, name):
        # create the dataset
         # first write boxlayout if it does not exist

        q = ldfab
        if len(q.FABs) > 0:
            ncomp = q.FABs[0].data.shape[-1::][0]
        else:
            ncomp = 0
        if self.IsOpen:
            if ncomp is 0:
                raise ValueError(" root has not fabs. Output should be assigned to a process with fabs")
            if 'boxes' not in self.Root[groupname].keys():
                BLO = [(id, (box.LoEnd() + box.HiEnd())) for id, box in zip(q.PIDs, q.Boxes)]
                self.WriteBoxLayout(BLO, groupname)

            if groupname not in self.Root.keys(): self.Root.create_group(groupname)
            if name not in self.Root[groupname].keys():

                size = 0
                offsets=list()
                for b in q.Boxes:
                    offsets.append(size)
                    size += ncomp * b.Ncells()

                offsets.append(size)


                # work on the attributes of the dataset
                self.Root[groupname].create_group(name + "_attributes")
                self.Root[groupname][name + "_attributes"].attrs['comps'] = ncomp
                b = [box for box in q.Boxes if q.FABs[0].box.ContainsBox(box)][0]
                ghosts=[s2>s1 for s1, s2 in zip(b.Size(), q.FABs[0].box.Size())]

                self.Root[groupname][name + "_attributes"].attrs['ghost'] = ghosts
                self.Root[groupname][name + "_attributes"].attrs['outputGhost'] = [0 for i in b.Size()]
                self.Root[groupname][name + "_attributes"].attrs['objectType'] = "b'FArrayBox'"

                #work on creating the datasets for value and offsets
                self.Root[groupname].create_dataset(name + ":datatype=0", (size,), compression="gzip")
                self.Root[groupname].create_dataset(name + ":offsets=0", data=offsets)


        for i,b in enumerate(q.Boxes):
            X=q.GatherBox(b, components=tuple(range(ncomp))).data.ravel(order='F')
            if self.IsOpen:
                begin = self.Root[groupname][name + ":offsets=0"][i]
                end = self.Root[groupname][name + ":offsets=0"][i+1]
                self.Root[groupname][name + ":datatype=0"][begin:end] = X





    # writes a fluxbox to the chkpoint file as '/groupname/'
    def WriteLevelDataFluxBox(self,ldvel, groupname, name):

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
                self.Root[groupname][name + "_attributes"].attrs['comps'] = ncomp
                b = [box for box in vel.Boxes if vel.FBs[0].box.ContainsBox(box)][0]
                ghosts=[s2>s1 for s1, s2 in zip(b.Size(), vel.FBs[0].box.Size())]

                self.Root[groupname][name + "_attributes"].attrs['ghost'] = ghosts
                self.Root[groupname][name + "_attributes"].attrs['outputGhost'] = [0 for i in b.Size()]
                self.Root[groupname][name + "_attributes"].attrs['objectType'] = "b'unknown'"

                #work on creating the datasets for value and offsets
                self.Root[groupname].create_dataset(name + ":datatype=0", (size,), compression="gzip")
                self.Root[groupname].create_dataset(name + ":offsets=0", data=offsets)


        for i,b in enumerate(vel.Boxes):
            X=vel.GatherFluxBox(b)
            if self.IsOpen:
                begin=self.Root[groupname][name+":offsets=0"][i]
                for fab in X.Fluxes:
                    S=1
                    for s in fab.data.shape:
                        S*=s
                    end=begin+S
                    self.Root[groupname][name+":datatype=0"][begin:end] = fab.data.ravel(order='F')
                    begin = end




    def WriteBoxLayout(self, BLO, group):
        if self.IsOpen:
            self.Root[group].create_dataset('boxes', data=[box for (id, box) in BLO])
            self.Root[group].create_dataset('Processors', data=[id for (id, box) in BLO])

class ReadCheckPoint:
    """ handles reading checkpoint files, providing SOMAR with info to buid LevelData<T> objects stored in the file"""
    def __init__(self, filename):
        self.FileHandle = h5.File(filename, 'r')
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
                    b = Box(box)
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
                    FBs.append(Args)

            Boxes = [b for i, b in self.BoxLayout[n]]
            ProcId=[(i,) for i,b in self.BoxLayout[n]]
            LevelDatas.append(LDFB((tuple(FBs), tuple(Boxes), tuple(ProcId))))

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
                    b = Box(box)
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
                    FABs.append(Args)

            Boxes = [b for i, b in self.BoxLayout[n]]
            ProcId=[(i,) for i,b in self.BoxLayout[n]]
            LDFABs = LD((tuple(FABs), tuple(Boxes), tuple(ProcId)))
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
            offsets = self.DataOffsets[n][FABName][:]
            Data = level[FABName + ":datatype=0"]

            # read the data
            for fluxbox in LevelDatas[n]:

                b = fluxbox.box

                for i, (id, valid) in enumerate(self.BoxLayout[n]):

                    ofst = offsets[i]
                    if b.ContainsBox(Box(valid)):
                        if id != comm.Get_rank():
                            raise ValueError("Something is wrong: this box should go to a different process")
                        for d in range(SpaceDim):
                            fab = fluxbox.Fluxes[d]
                            if SpaceDim == 3:
                                size = (fab.data.shape[0] -2) * (fab.data.shape[1] - 2) * (fab.data.shape[2] - 2)
                                S = (fab.data.shape[0]-2, fab.data.shape[1] - 2, fab.data.shape[2] - 2)
                                fab[1:-1, 1:-1, 1:-1, 0] = Data[ofst:ofst + size].reshape(S, order='F')
                                ofst+=size
                            else:
                                size = (fab.data.shape[0] -2) * (fab.data.shape[1] - 2)
                                S = (fab.data.shape[0]-2, fab.data.shape[1] - 2)
                                fab[1:-1, 1:-1, 0] = Data[ofst:ofst + size].reshape(S, order='F')
                                ofst += size
                        if ofst > offsets[i + 1]:
                            print(ofst, " ", offsets[i+1], " ", (offsets[i+1]-offsets[i]), " ", size , ' ', fab.data.shape )
                            raise ValueError("mismatch between data used and offsets.")




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
            if HasOutputGhosts:
                raise ValueError(FABName+ ": we cannot deal yet with outputghosts ")
            ncomps = self.LevelsGroupsAttributes[n][FABName+'_attributes']['comps']
            offsets = self.DataOffsets[n][FABName][:]
            Data = level[FABName + ":datatype=0"]

            # read the data
            for fab in LevelDatas[n]:
                b = fab.box

                for i, (id,valid) in enumerate(self.BoxLayout[n]):
                    if b.ContainsBox(Box(valid)):
                        if id != comm.Get_rank():
                            raise ValueError("Something is wrong: this box should go to a different process")
                        if HasOutputGhosts:
                            fab[...] = Data[offsets[i]:offsets[i + 1]].reshape(tuple(b.Size())+(ncomps,), order='F')
                        else:
                            if SpaceDim==3:
                                fab[1:-1, 1:-1, 1:-1,:] = Data[offsets[i]:offsets[i + 1]].reshape(tuple(Box(valid).Size()) + (ncomps,), order='F')
                            else:
                                fab[1:-1, 1:-1,:] = Data[offsets[i]:offsets[i + 1]].reshape(tuple(Box(valid).Size()) + (ncomps,), order='F')

    def IdAndBox(self, i, level):
        return self.BoxLayout[level][i]

    def NumberOfBoxes(self, level):
        return len(self.BoxLayout[level])





if __name__ == "__main__":
    dir3d=((1,0,0),(0,1,0),(0,0,1))
    R = ReadCheckPoint('chkpt_002902.2d.hdf5')
    level = 1
    level=min(level, len(R.LevelsAttributes[level]))
    print(R)
    qLevelDatas=R.SetUpLevelDatasFABs()
    R.ReadFABs(LevelDatas=qLevelDatas)
    pLevelDatas = R.SetUpLevelDatasFABs()
    R.ReadFABs(LevelDatas=pLevelDatas, FABName='pData')
    velLevelDatas=R.SetUpLevelDatasFluxBoxes()
    R.ReadFluxBoxes(LevelDatas=velLevelDatas)

    #Domain = Box((-128, -128, 0, 127, 127, 511))
    dom = R.LevelsAttributes[level]['domainBox']
    if len(dom) is 4:

        Domain = Box(dom)
    else:
        # just pick a slice midway in y
        y=(dom[4]+dom([1]))//2
        dom = (dom[0], y) + dom[2:3] + (y, dom[5])

    B = qLevelDatas[level].GatherBox(Domain, components=(3,))
    if (MPI.COMM_WORLD.Get_rank() == 0):
        plt.imshow(B[:,:,0].transpose())
        plt.colorbar()
        plt.savefig("B.pdf")
        plt.close()
    V=velLevelDatas[level].GatherFluxBox(Domain, FillWith=0.0)

    if (MPI.COMM_WORLD.Get_rank() == 0):
        for fluxes in V.Fluxes:
            plt.imshow(fluxes[:,:,0].transpose())
            plt.colorbar()
            plt.savefig("U.pdf")
            plt.close()
            plt.imshow(fluxes[:,:,0].transpose())
            plt.colorbar()
            plt.savefig("V.pdf")
            plt.close()


    # velLevelDatas=R.SetUpLevelDatasFluxBoxes()
    # R.ReadFluxBoxes(LevelDatas=velLevelDatas)
    #if MPI.COMM_WORLD.Get_rank()==0: print(R)
    # for fun, calculate vorticity
    dx=R.LevelsAttributes[0]['vec_dx']
    divmax=0.
    size=[s for s in V.box.Size()]
    div=np.zeros(size,dtype='double', order='F')
    if MPI.COMM_WORLD.Get_rank()==0:
       # print(size)

        for d,fab in enumerate(V.Fluxes):

            div[1:-1,1:-1]+=(fab[1:-1,1:-1,0]-fab[1-dir3d[d][0]:-dir3d[d][0]-1,
                                                            1 - dir3d[d][1]: - dir3d[d][1] - 1,
                                                            0])/dx[d]



        divmax=max(divmax,np.max(np.abs(div)))
        print("maximum divergence ", divmax)
    # del div
    # del V
    # del B
    # if MPI.COMM_WORLD.Get_rank() ==0 :
    #     print(R)

    # if MPI.COMM_WORLD.Get_rank() == 0:
    #     f = h5.File('testcheck.hdf5', 'w')

    #     Groups = [f.create_group(name) for name in R.GroupNames]
    #     DomainBox=[Box(LevelAttributes['domainBox']) for LevelAttributes in R.LevelsAttributes]
    # # work on the root attributes.
    #     for (name, value) in R.RootAttributes.items():
    #         f['/'].attrs.create(name, value)
    #     # work on group attributes
    #     for (name, value) in R.ChomboGlobalAttributes.items():
    #         Groups[0].attrs.create(name, value)
    #     # work on level attributes
    #     for i, level in enumerate(R.Levels):
    #         for (name, value) in R.LevelsAttributes[i].items():
    #             Groups[i + 1].attrs.create(name, value)

    #     # work on boxlayout

    #     for i, level in enumerate(R.Levels):
    #         Groups[i + 1].create_dataset('boxes', data=[box for (id, box) in R.BoxLayout[i]])
    #         Groups[i + 1].create_dataset('Processors', data=[id for (id, box) in R.BoxLayout[i]])
    #         [Groups[i + 1].create_dataset(name+':offsets=0', data=value) for (name, value) in R.DataOffsets[i].items()]




    #     # work on level subgroups and datasets
    #     for i, levelGroupsAttributes in enumerate(R.LevelsGroupsAttributes):
    #         for key in levelGroupsAttributes.keys():
    #             Groups[i + 1].create_group(key)
    #             [Groups[i+1][key].attrs.create(name, value) for (name,value) in levelGroupsAttributes[key].items()]

    #         for name,ofst in R.DataOffsets[i].items():
    #             size = (ofst[-1::][0],)
    #             # size = [R.LevelsGroupsAttributes[i][name+"_attributes"]['comps']]+DomainBox[i].Size()[::-1]
    #             # print("size of ", name, " is ",size)

    #             Groups[i+1].create_dataset(name + ':datatype=0', size, compression="gzip")

    #     # collect the boxes and write them to the file
    # for n in range(len(R.Levels)):
    #     for i,(id, box) in enumerate(R.BoxLayout[n]):
    #         valid = Box(box)


    #         ncomp=R.LevelsGroupsAttributes[n]['qData_attributes']['comps']
    #         X = qLevelDatas[n].GatherBox(valid, components=tuple(range(ncomp))).data.ravel(order='F')

    #         if MPI.COMM_WORLD.Get_rank() == 0:

    #             begin = R.DataOffsets[n]['qData'][i]
    #             end = R.DataOffsets[n]['qData'][i + 1]
    #             Groups[n + 1]['qData:datatype=0'][begin:end] = X


    #             #note the C order below!
    #             # begin = [s - c for s, c in zip(valid.LoEnd(), DomainBox[n].LoEnd())][::-1]
    #             # end = [s - c + 1 for s, c in zip(valid.HiEnd(), DomainBox[n].LoEnd())][::-1]
    #             # print(i,begin, end)
    #             # Groups[n + 1]['qData:datatype=0'][:,
    #             #                                   begin[2]: end[2],
    #             #                                   begin[1]: end[1],
    #             #                                   begin[0]: end[0],
    #             #                                   ] = X[...].ravel(order='F').reshape((8,64,64,64), order='C')
    #         ncomp=R.LevelsGroupsAttributes[n]['pData_attributes']['comps']
    #         X = pLevelDatas[n].GatherBox(valid, components=tuple(range(ncomp))).data.ravel(order='F')
    #         if MPI.COMM_WORLD.Get_rank() == 0:
    #             begin = R.DataOffsets[n]['pData'][i]
    #             end=R.DataOffsets[n]['pData'][i+1]
    #             Groups[n + 1]['pData:datatype=0'][begin:end] = X
    #             # begin = [s - c for s, c in zip(valid.LoEnd(), DomainBox[n].LoEnd())][::-1]
    #             # end = [s - c + 1 for s, c in zip(valid.HiEnd(), DomainBox[n].LoEnd())][::-1]
    #             # Groups[n + 1]['pData:datatype=0'][:,
    #             #                                   begin[2]: end[2],
    #             #                                   begin[1]: end[1],
    #             #                                   begin[0]: end[0],
    #             #                                   ] = X[...].ravel(order='F').reshape((1,64,64,64), order='C')
    #         X=velLevelDatas[n].GatherFluxBox(valid)
    #         if MPI.COMM_WORLD.Get_rank( ) == 0:
    #             begin=R.DataOffsets[n]['velData'][i]
    #             for fab in X.Fluxes:
    #                 S=1
    #                 for s in fab.data.shape:
    #                     S*=s
    #                 end=begin+S
    #                 Groups[n + 1]['velData:datatype=0'][begin:end] = fab.data.ravel(order='F')
    #                 begin=end


    #         # close the file
    # if MPI.COMM_WORLD.Get_rank()==0:
    #     f.close()
    # MPI.COMM_WORLD.barrier()
    # S = ReadCheckPoint('testcheck.hdf5')
    # if MPI.COMM_WORLD.Get_rank() == 0:
    #     print(S)


    # LDs=R.SetUpLevelDatasFABs(FABName='qData')
    # R.ReadFABs(LevelDatas=LDs, FABName='qData')

    # LDFBs = R.SetUpLevelDatasFluxBoxes()
    # R.ReadFluxBoxes(LevelDatas=LDFBs)


    # Domain = Box(LDs[0].PhysDom()[0] + LDs[0].PhysDom()[1])

    # Domain=Box((-128,0,0,127,0,511))
    # # print(Domain)

    # X = LDs[0].GatherBox(Domain, components=tuple(range(8)))

    # if MPI.COMM_WORLD.Get_rank() == 0:
    #     for c in range(X.data.shape[3]):
    #         plt.imshow(X[:,0, :,c].transpose())
    #         plt.colorbar()
    #         plt.savefig(str(c)+'test.pdf')
    #         plt.close()



