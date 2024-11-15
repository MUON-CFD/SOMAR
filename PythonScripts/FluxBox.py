import PyGlue as pg
import numpy as np
from Box import Box
from FAB import FAB
class FluxBox:
    """ define a FluxBox as a tuple of FABs, plus metadata of the common box (box, centering, number of components)"""

    def __init__(self, Args, FillWith=None):
        """construct the FluxBox specified by Args       : Args[0]-> tuple of lower index of box
                                                       Args[1]-> tuple of size
                                                       Args[2]-> tuple of centering
                                                       Args[3]-> number of components
                                                       Args[4] -> FAB of first comp
                                                       Args[5] -> FAB of second comp... """
        if Args[-1] != "FluxBox":
            raise ValueError("FluxBox constructor called with non FluxBox Args")

        self.size = Args[1]
        self.dimension = len(self.size)

        self.nComp = Args[3][0]
        if Args[2] is not None:
            self.centering = Args[2]
        else:
            self.centering = tuple([0 for x in self.size])
        dummy = list(Args[0])
        for (l,s) in zip(Args[0], self.size):
            dummy.append(l + s - 1)
        dummy.append('Box')
        self.box = Box(dummy)
        Fluxes=[]
        if len(Args[:-1]) > 4 and all([x is not None for x in Args[4:-1]]):
            for fab in Args[4:-1]:

                Fluxes.append(FAB(fab))

        else:

            for d in range(self.dimension):
                fab=[]
                dummy = [s  for s in Args[0]]
                #dummy[d] += 1
                fab.append(dummy)
                dummy = [s  for s in Args[1]]
                dummy[d] += 1

                fab.append(dummy)
                dummy = [s for s in Args[2]]
                dummy[d] = 1
                fab.append(dummy)
                fab.append((self.nComp,))
                fab.append(None)
                fab.append('FAB')
                Fluxes.append(FAB(fab, FillWith=FillWith))

        self.Fluxes = tuple(Fluxes)

    def __iter__(self):
        return iter(self.Fluxes)

    def __next__(self):
        return next(self.Fluxes)

    def __str__(self):
        axis={0:'x', 1:'y', 2:'z'}
        string=""
        for d in range(self.dimension):
                string += "Flux along " + axis[d] + "-direction \n"
                string += str(self.Fluxes[d]) + "\n"

        return string

    def copyFrom(self, srcFB, box=None, destComp=(0,), srcComp=(0,)):
        IB=self.box.Intersection(srcFB.box)
        if box is not None:
            IB = IB.Intersection(box)

        for d, srcFab in enumerate(srcFB.Fluxes):
            LoEnd=list(IB.LoEnd())
            Size=list(IB.Size())
            LoEnd[d]+=1
            Size[d]+=0
            B=Box.FromLoEndAndSize(LoEnd,Size)

            self.Fluxes[d].copyFrom(srcFab, box=B, destComp=destComp, srcComp=srcComp)


    def Write(self,filename='FB.npz'):
        if self.dimension==2:
            np.savez(filename, U=self.Fluxes[0].data, V=self.Fluxes[1].data)
        else:
            np.savez(filename, U=self.Fluxes[0].data, V=self.Fluxes[1].data, W=self.Fluxes[2].data)

    def Read(self, filename='FB.npz'):
        with np.load(filename) as data:
            try:
                self.Fluxes[2].data=data['W']
            except:
                pass
            finally:
                self.Fluxes[0].data=data['U']
                self.Fluxes[1].data=data['V']









