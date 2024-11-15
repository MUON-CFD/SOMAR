import numpy as np
import PyGlue as pg
from Box import Box as Box
from FAB import FAB as FAB
from mpi4py import MPI

class LevelDataFAB:



    """ Handles LevelData<FArrayBox> """
    def __init__(self, Arg, offsets=None):
        self.comm = MPI.COMM_WORLD
        self.MpiNumProcs = self.comm.Get_size()
        self.myMpiId = self.comm.Get_rank()
        self.tags=[]

        self.FABs, self.Boxes, self.PIDs, self.Domain = self._LDtoNumpy(Arg, offsets)

        self.SizesForIO, self.tags=self._SizesForIO()
        self.SpaceDim=len(self.Boxes[0].LoEnd())
        try:
            nComp=self.FABs[0].data.shape[-1]
        except:
            nComp=None
        finally:
            self.nComp=MPI.COMM_WORLD.bcast(nComp, root=0)

    @classmethod 
    def MakeFromLDFAB(cls, ldfab):
        """ returns an LDFAB defined on the same grids as ldfab. New memory is allocated"""
        Args = 5*[None]
        Args[4] = 'LevelDataFAB'
        Args[3] = tuple(ldfab.Domain.LoEnd()+ldfab.Domain.HiEnd()+('Box',))
        Args[2] = [(id,) for id in ldfab.PIDs]
        Args[1] = [(b.LoEnd() + b.HiEnd()+('Box',)) for b in ldfab.Boxes]
        Args[0] = [FAB.FromFAB(fab).FAB2C() for fab in ldfab.FABs]
        return cls(Args)

    #@classmethod
    #def MakeFromF
    def __iter__(self):
        """ the iterator of the class iterates over the FABs"""
        return iter(self.FABs)

    def __next__(self):
        return next(self.FABs)

    def __str__(self):
        string = str(self.SpaceDim)+"-dimensional LevelDataFAB \n"
        string += "=================================\n"
        string += "My rank is " + str(self.myMpiId) + " out of a total of " + str(self.MpiNumProcs) + "\n"
        string += "=================================\n"
        string += " Box Layout \n"
        for Id,box in zip(self.PIDs,self.Boxes):
            string += " Proc Id " + str(Id) + " " + str(box) + "\n"
        string += "=================================\n"
        string += "Local fields \n"
        for i, fab in enumerate(self):
            string += " field # " + str(i) + "\n"
            string += str(fab)
        return string

    def _SizesForIO(self):
        # calculate space requirement (to be used in IO)
        Sizes=dict()
        tags=[]
        for q in self:
            try:
                box=[box for box in self.Boxes if q.box.ContainsBox(box)][0]
            except:
                raise ValueError("we have a fluxBox without a box in the list of boxes")


            tag=str(hash(box.LoEnd() ))
            Sizes[tag]=q[...].shape
            tags.append(tag)
        # combine in a single dictionary
        #return self.comm.allgather(Sizes)
        return dict((key,value[key]) for value in self.comm.allgather(Sizes) for key in value), tuple(tags)

    def LevelDataFAB2C(self):
        """ returns a tuple that can be used to initialize a new LevelData which aliases this one"""
        Args = 5*[None]
        Args[4] = 'LevelDataFAB'
        Args[3] = tuple(self.Domain.LoEnd()+self.Domain.HiEnd()+('Box',))
        Args[2] = [(id,) for id in self.PIDs]
        Args[1] = [(b.LoEnd() + b.HiEnd()+('Box',)) for b in self.Boxes]
        Args[0] = [fab.FAB2C() for fab in self.FABs]
        return Args
    

    def PhysDom(self):
        """ Determines the extent of the domain from a leveldata. Note that this may not the total extent of the
            actual domain. It is just the smallest box that contains all the boxes in the LD collection """
        lo = [x for x in self.Boxes[0].LoEnd()]
        hi = [x for x in self.Boxes[0].HiEnd()]
        for Box in self.Boxes:
            for d in range(len(lo)):
                lo[d] = min(lo[d], Box.LoEnd()[d])
                hi[d] = max(hi[d], Box.HiEnd()[d])
        return (tuple(lo), tuple(hi))

    def GatherBox(self, box, components=(0,), RootId=0, FillWith=np.nan):
        """ Given a box and a tuple of components, fills a FAB on myMpiId==0 with data taken from this LevelDataFAB. Only
        the FAB on RootId is filled with real data. Everybody else returns a FAB filled with nans"""

        Args = []
        Args.append(box.LoEnd())
        Args.append(box.Size())
        Args.append(box.Centering())
        Args.append((len(components),))
        Args.append(None)  # we want Python to allocate memory
        Args.append('FAB')
        X = FAB(Args, FillWith=FillWith)

        # list of boxes in the LevelData that have a nonzero intersection with box
        Active = [(Box, PID) for (Box, PID) in zip(self.Boxes, self.PIDs) if not Box.Intersection(box).IsEmpty()]

        # of those, this is a list that the process owns
        MyActive = [Box for (Box, pid) in Active if self.myMpiId == pid]

        for nc,c in enumerate(components):
            # nonroot processes send
            if self.myMpiId != RootId:
                # the double for loops over the intersection of FABs on this process that intersect with box
                for fab in self:
                    if any([fab.box.ContainsBox(b1) for b1 in MyActive]):
                        b2,_=[ (Box,PID) for (Box, PID) in Active if fab.box.ContainsBox(Box)][0]
                        intBox = fab.box.Intersection(box)
                        intBox = intBox.Intersection(b2)
                        if (intBox.IsEmpty()):
                            raise ValueError(" Something is wrong, we are having an empty intersection ")

                        Args = []
                        Args.append(intBox.LoEnd())
                        Args.append(intBox.Size())
                        Args.append(intBox.Centering())
                        Args.append((1 ,))
                        Args.append(None)  # we want Python to allocate memory
                        Args.append('FAB')
                        dummyFAB = FAB(Args, FillWith=np.nan)
                        dummyFAB.copyFrom(fab, srcComp=(c,))

                        # generate unique tag
                        tag=hash(dummyFAB.box.LoEnd())%10000
                        req=self.comm.Isend(dummyFAB.data, dest=RootId, tag=tag) # use non blocking sends!
                        req.Wait()

            else:  # process where gathers
                # first from himself
                # the double for loops over the intersection of FABs on this process that intersect with box
                for fab in self:
                    if any([fab.box.ContainsBox(b1) for b1 in MyActive]):
                        b2=[Box for (Box, PID) in Active if fab.box.ContainsBox(Box)][0]
                        intBox = fab.box.Intersection(box)
                        intBox = intBox.Intersection(b2)

                        if (intBox.IsEmpty()):
                            raise ValueError(" Something is wrong, we are having an empty intersection ")

                        X.copyFrom(fab, box=intBox, destComp=(nc,), srcComp=(c,))

                for (Box, SenderId) in Active:
                    if SenderId != RootId:
                        intBox = Box.Intersection(box)
                        Args = []
                        Args.append(intBox.LoEnd())
                        Args.append(intBox.Size())
                        Args.append(intBox.Centering())
                        Args.append((1 ,))
                        Args.append(None)  # we want Python to allocate memory
                        Args.append('FAB')
                        dummyFAB = FAB(Args, FillWith=np.nan)
                        tag=hash(dummyFAB.box.LoEnd())%10000
                        self.comm.Recv(dummyFAB.data, source=SenderId, tag=tag)

                        X.copyFrom(dummyFAB, destComp=(nc,), srcComp=(0,))
        return X

    def _LDtoNumpy(self, Args, offsets):
        """ unpacks a LevelData into a tuple of tuples"""
        FABs = []
        if Args[-1] != "LevelDataFAB":
            raise ValueError("LevelDataFAB called with non LevelDataFAB ARgs")

        for F in Args[0]:

            FABs.append(FAB(F,offsets=offsets))

        Boxes = []
        for F in Args[1]:

            Boxes.append(Box(F))

        PIDs=tuple([p[0] for p in Args[2]])

        Domain = Box(Args[3])
        
        return (tuple(FABs), tuple(Boxes), PIDs, Domain)


if __name__ == "__main__":
    import Slicer
    from random import seed
    from random import randint
    import matplotlib.pyplot as plt
    # seed random number generator
    seed(1)
    count=0
    nx, ny, nz = (256, 128, 512) # total domain
    mx, my, mz = (256//4, 128//8, 512//32)  # min block
    while(count<1):
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

        LD = LevelDataFAB(Args)

        #fill the fabs with unique markers

        for f in LD:
            f[...] = -10 * (MyRank + 1)
            b=f.box
            start=(b.LoEnd()[0]+1)*mx+(b.LoEnd()[1]+1)*my+(b.LoEnd()[2]+1)*mz
            f[1:-1, 1:-1, 1:-1, 0] = np.arange(start,start+mx * my * mz,1).reshape((mx, my, mz), order='F')






        B = Box.FromLoEndAndSize((0, 0, 10), (nx, ny, 1))
        # work on an aliased copy
        LDAlias = LevelDataFAB(LD.LevelDataFAB2C())

        X = LDAlias.GatherBox(B)
        Y=Slicer.GatherSlice(LDAlias, dir=2, pos=10)
        if MyRank == 0:
            #print(X[126:135,0])

            # check consistency
            Active = [(b,id) for b, id in zip(LD.Boxes, LD.PIDs) if not b.Intersection(B).IsEmpty()]


            for (b, id) in Active:

                Args = []
                Args.append(tuple(b.LoEnd()))
                Args.append(tuple(b.Size()))
                Args.append((0, 0, 0))
                Args.append((1,))
                Args.append(None)
                Args.append('FAB')
                f = FAB(Args)
                b=f.box
                start=b.LoEnd()[0]*mx+b.LoEnd()[1]*my+b.LoEnd()[2]*mz
                f[:,:,:,0] = np.arange(start,start+mx * my * mz,1).reshape((mx, my, mz), order='F')

                f[...] *= -1



                X.Plus(f)
                Y.Plus(f)

            print(" iteration ", count, "Min err ", np.min(X.data), "and Max err ", np.max(X.data))
            print(" iteration ", count, "Min err slicer", np.min(Y.data), "and Max err ", np.max(Y.data))
            # print(np.max(np.abs(X.data)))
            # print(X.data.shape)
            # print(X[126:135,0])
            dir=2
            pos=10
            filename='Q_iter_'+str(count).zfill(7)+'_dir_'+str(dir)+'_pos_'+str(pos)+'.npy'

            Y.Write()


            X.data[...]=1e10
            X.Read()

            print(X[1,2,0], Y[1,2,0])
        count += 1
            # if np.max(np.abs(X.data)) > 0:

            #     plt.imshow(X.data[:,:, 0, 0].transpose())
            #     plt.colorbar()
            #     plt.show()







