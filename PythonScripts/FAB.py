import PyGlue as pg
import numpy as np
from Box import Box
class FAB:
    """ define a FAB, encapsulating data and metadata (box, centering, number of components). It also adds defines some
        useful operations. In particular, FAB[key] accesses FAB.data[key]"""

    def __init__(self, Args, FillWith=None, offsets=None):
        """construct the FAB specified by Args       : Args[0]-> tuple of lower index of box
                                                       Args[1]-> tuple of size
                                                       Args[2]-> tuple of centering
                                                       Args[3]-> number of components
                                                       Args[4] -> address of buffer. If None, a new numpy is allocated."""
        if Args[-1] != "FAB":
            raise ValueError("FAB constructor called on a non FAB argument")
        self.size = Args[1]

        self.ncomp = Args[3][0]
        if Args[2] is not None:
            self.centering = Args[2]
        else:
            self.centering = tuple([0 for x in self.size])
        dummy = list(Args[0])
        for (l,s) in zip(Args[0], self.size):
            dummy.append(l + s - 1)
        dummy.append('Box')
        self.box = Box(dummy)
        if Args[4] is not None:
            self.data = pg.FABToNumpy(Args)
        else:
            shape = list(self.size)
            shape.append(self.ncomp)

            shape=tuple(shape)
            self.data = np.empty(shape, dtype='double', order='F')
            if FillWith is not None:
                self.data[...] = FillWith
        self.dimension = len(self.size)
        if offsets is not None:
            self.offsets = offsets
        else:
            self.offsets = self.box.LoEnd()

    def __str__(self):
        return 'FAB with '+str(self.ncomp)+' components defined on '+str(self.box)

    def __getitem__(self,key):
        return self.data[key]

    def __setitem__(self,key,value):
        self.data[key]=value

    def FAB2C(self):
        """ returns a tuple that can be used to create a new FAB which aliases the numpy of this one. """
        Args = 6 * [None]
        Args[0] = self.box.LoEnd()
        Args[1] = self.box.Size()
        Args[2] = self.centering
        Args[3] = (self.ncomp,)
        Args[4] = self.data.view()
        Args[5] = 'FAB'
        return Args

    @classmethod
    def FromBox(cls, box, ncomp=1):
        Args = []
        Args.append(box.LoEnd())
        Args.append(box.Size())
        Args.append(box.Centering())
        Args.append((ncomp ,))
        Args.append(None)  # we want Python to allocate memory
        Args.append('FAB')
        return cls(Args, FillWith=np.nan)
    
    @classmethod
    def FromFAB(cls, FAB):
        Args=FAB.FAB2C()
        Args[4]=None
        return cls(Args)

    def Plus(self, srcFAB, destComp=(0,), srcComp=(0,)):
        """ adds srcFAB to self on the intersection of their boxes.
        The components of the source specified by srcComp are added
        to the destComp. """

        for (dC, sC) in zip(destComp, srcComp):
            IB = self.box.Intersection(srcFAB.box)

            if self.dimension == 3:
                self.data[IB.LoEnd()[0] - self.offsets[0]: IB.HiEnd()[0] - self.offsets[0] + 1,
                          IB.LoEnd()[1] - self.offsets[1]: IB.HiEnd()[1] - self.offsets[1] + 1,
                          IB.LoEnd()[2] - self.offsets[2]: IB.HiEnd()[2] - self.offsets[2] + 1, dC] += \
                srcFAB.data[IB.LoEnd()[0] - srcFAB.offsets[0]: IB.HiEnd()[0] - srcFAB.offsets[0] + 1,
                            IB.LoEnd()[1] - srcFAB.offsets[1]: IB.HiEnd()[1] - srcFAB.offsets[1] + 1,
                            IB.LoEnd()[2] - srcFAB.offsets[2]: IB.HiEnd()[2] - srcFAB.offsets[2] + 1, sC]

            else:
                self.data[IB.LoEnd()[0] - self.offsets[0]: IB.HiEnd()[0] - self.offsets[0] + 1,
                          IB.LoEnd()[1] - self.offsets[1]: IB.HiEnd()[1] - self.offsets[1] + 1, dC] += \
                srcFAB.data[IB.LoEnd()[0] - srcFAB.offsets[0]: IB.HiEnd()[0] - srcFAB.offsets[0] + 1,
                            IB.LoEnd()[1] - srcFAB.offsets[1]: IB.HiEnd()[1] - srcFAB.offsets[1] + 1, sC]

    def copyFrom(self, srcFAB, box=None, destComp=(0,), srcComp=(0,)):
        """ copies srcFAB to self on the intersection of their boxes with box (if present).
        The components of the source specified by srcComp are added
        to the destComp. """


        for (dC, sC) in zip(destComp, srcComp):

            IB = self.box.Intersection(srcFAB.box)

            if box is not None:
                IB = IB.Intersection(box)

            if self.dimension == 3:

                self.data[IB.LoEnd()[0] - self.offsets[0]: IB.HiEnd()[0] - self.offsets[0] + 1,
                          IB.LoEnd()[1] - self.offsets[1]: IB.HiEnd()[1] - self.offsets[1] + 1,
                          IB.LoEnd()[2] - self.offsets[2]: IB.HiEnd()[2] - self.offsets[2] + 1, dC] = \
                srcFAB.data[IB.LoEnd()[0] - srcFAB.offsets[0]: IB.HiEnd()[0] - srcFAB.offsets[0] + 1,
                            IB.LoEnd()[1] - srcFAB.offsets[1]: IB.HiEnd()[1] - srcFAB.offsets[1] + 1,
                            IB.LoEnd()[2] - srcFAB.offsets[2]: IB.HiEnd()[2] - srcFAB.offsets[2] + 1, sC]

            else:
                self.data[IB.LoEnd()[0] - self.offsets[0]: IB.HiEnd()[0] - self.offsets[0] + 1,
                          IB.LoEnd()[1] - self.offsets[1]: IB.HiEnd()[1] - self.offsets[1] + 1, dC] = \
                srcFAB.data[IB.LoEnd()[0] - srcFAB.offsets[0]: IB.HiEnd()[0] - srcFAB.offsets[0] + 1,
                            IB.LoEnd()[1] - srcFAB.offsets[1]: IB.HiEnd()[1] - srcFAB.offsets[1] + 1, sC]


    def Write(self,filename='FAB.npy'):
        np.save(filename, self.data)



    def Read(self,filename='FAB.npy'):
        self.data=np.load(filename)

    def retrieve(self):
        pass


