import numpy as np
import PyGlue as pg
from Box import Box
from FluxBox import FluxBox
from FAB import FAB
from mpi4py import MPI
class LevelDataFluxBox:



    """ Handles LevelData<FluxBox> """
    def __init__(self, Arg, offsets=None):
        self.comm = MPI.COMM_WORLD
        self.MpiNumProcs = self.comm.Get_size()
        self.myMpiId = self.comm.Get_rank()
        self.FBs, self.Boxes, self.PIDs , self.Domain= self._LDtoNumpy(Arg, offsets)
        self.SpaceDim = len(self.Boxes[0].LoEnd())
        try:
            nComp=self.FBs[0].nComp
        except:
            nComp=None
        finally:
            self.nComp=self.comm.bcast(nComp, root=self.PIDs[0])
        self.SizesForIO,self.tags=self._SizesForIO()


    def __iter__(self):
        """ The iterator of the class iterates over the FluxBoxes"""
        return iter(self.FBs)

    def __next__(self):
        return next(self.FBs)

    def __str__(self):
        string = str(self.SpaceDim)+"-dimensional LevelDataFluxBox \n"
        string += "=================================\n"
        string += "My rank is " + str(self.myMpiId) + " out of a total of " + str(self.MpiNumProcs) + "\n"
        string += "=================================\n"
        string += " Box Layout \n"
        for Id,box in zip(self.PIDs,self.Boxes):
            string += " Proc Id " + str(Id) + " " + str(box) + "\n"
        string += "=================================\n"
        string += "Local fields \n"
        for i, flux in enumerate(self):
            string += " field # " + str(i) + "\n"
            string += str(flux)
        return string

    def _SizesForIO(self):
        # calculate space requirement (to be used in IO)
        Sizes=dict()
        tags=list()
        for q in self:
            try:
                box=[box for box in self.Boxes if q.box.ContainsBox(box)][0]
            except:
                raise ValueError("we have a fluxBox without a box in the list of boxes")
            tag=self.SpaceDim*[None]
            for d,flux in enumerate(q.Fluxes):
                tag[d]=str(hash(box.LoEnd()+(d,)) )
                Sizes[tag[d]]=flux[...].shape
            tags.append(tuple(tag))
        # combine in a single dictionary
        return dict((key,value[key]) for value in self.comm.allgather(Sizes) for key in value), tuple(tags)

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
    def GatherFluxBox(self,box, RootId=0, FillWith=np.nan):
        """ given a box, fills a FluxBox on RootIt with a fluxbox on the box. Only RootId returns X filled with valid
        data. All other processes get NaN"""

        Args=[]
        Args.append(box.LoEnd())
        Args.append(box.Size())
        Args.append(box.Centering())

        Args.append((1,))
        for d in range(len(box.LoEnd())): Args.append(None)
        Args.append('FluxBox')
        X=FluxBox(Args,FillWith=FillWith)

        components=(0,)
        # list of boxes in the LevelData that have a nonzero intersection with box

        Active = [(Box, PID) for (Box, PID) in zip(self.Boxes, self.PIDs) if not Box.Intersection(box).IsEmpty()]

        # of those, this is a list that the process owns
        MyActive = [Box for (Box, pid) in Active if self.myMpiId == pid]

        for nc,c in enumerate(components):
            # nonroot processes send
            if self.myMpiId != RootId:
                # the double for loops over the intersection of FABs on this process that intersect with box
                for fb in self:
                    if any([fb.box.ContainsBox(b1) for b1 in MyActive]):
                        b2,_=[ (Box,PID) for (Box, PID) in Active if fb.box.ContainsBox(Box)][0]
                        intBox = fb.box.Intersection(box)
                        intBox = intBox.Intersection(b2)
                        if (intBox.IsEmpty()):
                            raise ValueError(" Something is wrong, we are having an empty intersection ")
                        for d in range(len(intBox.LoEnd())):
                            Args = []
                            Args.append(list(intBox.LoEnd()))
                            Args[0][d]+=1
                            Args.append(list(intBox.Size()))
                            Args[1][d]+=0
                            Args.append(intBox.Centering())
                            Args.append((self.nComp ,))
                            Args.append(None)  # we want Python to allocate memory
                            Args.append('FAB')
                            dummyFAB = FAB(Args, FillWith=np.nan)
                            dummyFAB.copyFrom(fb.Fluxes[d], srcComp=(c,))


                            # generate unique tag

                            tag=hash(dummyFAB.box.LoEnd())%10000
                            req=self.comm.Isend(dummyFAB.data, dest=RootId, tag=tag) # use non blocking sends!
                            req.Wait()

            else:  # process where gathers
                # first from himself
                # the double for loops over the intersection of FABs on this process that intersect with box
                for fb in self:
                    if any([fb.box.ContainsBox(b1) for b1 in MyActive]):
                        b2=[Box for (Box, PID) in Active if fb.box.ContainsBox(Box)][0]
                        intBox = fb.box.Intersection(box)
                        intBox = intBox.Intersection(b2)
                        if (intBox.IsEmpty()):
                            raise ValueError(" Something is wrong, we are having an empty intersection ")

                        X.copyFrom(fb, box=intBox, destComp=(nc,), srcComp=(c,))

                for (Box, SenderId) in Active:
                    if SenderId != RootId:
                        intBox = Box.Intersection(box)
                        FBArgs=list()
                        FBArgs.append(intBox.LoEnd())
                        FBArgs.append(intBox.Size())
                        FBArgs.append(intBox.Centering())
                        FBArgs.append((1,))

                        for d in range(len(intBox.LoEnd())):
                            Args = []
                            Args.append(list(intBox.LoEnd()))
                            Args[0][d]+=1
                            Args.append(list(intBox.Size()))
                            Args[1][d]+=0
                            Args.append(intBox.Centering())
                            Args.append((1 ,))
                            Args.append(None)  # we want Python to allocate memory
                            Args.append('FAB')
                            dummyFAB = FAB(Args, FillWith=np.nan)
                            tag=hash(dummyFAB.box.LoEnd())%10000
                            self.comm.Recv(dummyFAB.data, source=SenderId, tag=tag)
                            Args[4]=dummyFAB.data.view()
                            FBArgs.append(Args)
                        FBArgs.append('FluxBox')
                        dummyFluxBox=FluxBox(FBArgs)

                        X.copyFrom(dummyFluxBox, destComp=(nc,), srcComp=(0,))

        return X

    def _LDtoNumpy(self, Args, offsets):
        """ unpacks a LevelData into a tuple of tuples"""
        if Args[-1] != "LevelDataFluxBox":
            raise ValueError("LevelDataFluxBox constructed with non LevelDataFluxBox Args")
        FBs = []

        for F in Args[0]:
            FBs.append(FluxBox(F))

        Boxes = []
        for F in Args[1]:

            Boxes.append(Box(F))

        PIDs = tuple([p[0] for p in Args[2]])
       
        Domain = Box(Args[3])
        return (tuple(FBs), tuple(Boxes), PIDs, Domain)


