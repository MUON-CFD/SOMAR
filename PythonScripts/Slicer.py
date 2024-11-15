

import numpy as np
from LevelDataFAB import LevelDataFAB
from LevelDataFB import LevelDataFluxBox
from Box import Box


"""
In 2D Slice and Line both return 1D objects... 
"""
def GatherSlice(LD,dir=0, pos=0):
    """
    Returns a FAB or FluxBox (depending on the input) sliced over a plane with normal
    in direction dir and at logical position given by pos. If the plane is outside the
    enveloping box, returns None. Note that when the plane is within the enveloping box,
    only the root process will have actual data.
    """
    if dir>LD.SpaceDim-1:
        raise ValueError("dir must be less than SpaceDim-1")

    LoEnd,HiEnd=LD.PhysDom()
    DomBox=Box(LoEnd+HiEnd+('Box',))
    LoEnd=list(LoEnd)
    HiEnd=list(HiEnd)
    LoEnd[dir]=pos
    HiEnd[dir]=pos
    SliceBox=Box(tuple(LoEnd)+tuple(HiEnd)+('Box',))
    if DomBox.ContainsBox(SliceBox):
    # note that only the proc with procID=0 will have data. All other processes
    # have nans.
        try:
            return LD.GatherBox(SliceBox,components=[i for i in range(LD.nComp)])
        except:
            return LD.GatherFluxBox(SliceBox)
    else:
        return None

def GatherLine(LD, dir=0, pos=(0,0,0)):
    """
     Returns a FAB or FluxBox (depending on the input) sliced over a line
    in direction dir passing by the logical  pos in the normal plane. Note that only the values of pos
    not in direction dir are used. If the line is outside the
    enveloping box, returns None. Note that when the line is within the enveloping box,
    only the root process will have actual data.
    """

    
    if dir>LD.SpaceDim-1:
        raise ValueError("dir must be less than SpaceDim-1")
    if len(pos)!=LD.SpaceDim:
        raise ValueError("len(pos) must be SpaceDim")
    LoEnd,HiEnd=LD.PhysDom()
    DomBox=Box(LoEnd+HiEnd+('Box',))
    LoEnd=list(LoEnd)
    HiEnd=list(HiEnd)
    for d in range(LD.SpaceDim):
        if d!=dir:
            LoEnd[d]=pos[d]
            HiEnd[d]=pos[d]

    SliceBox=Box(tuple(LoEnd)+tuple(HiEnd)+('Box',))

    if DomBox.ContainsBox(SliceBox):
        
    # note that only the proc with procID=0 will have data. All other processes
    # have nans.
        try:
        
            return LD.GatherBox(SliceBox,components=[i for i in range(LD.nComp)])
        
        except:
        
            return LD.GatherFluxBox(SliceBox)
    else:
        print("we should not be here")
        return None

def GatherPoint(LD, pos):
    """
     Returns a FAB or FluxBox (depending on the input) sliced over a point at position pos
     If the point is outside the
    enveloping box, returns None. Note that when the line is within the enveloping box,
    only the root process will have actual data.
    """
    if len(pos)!=LD.SpaceDim:
        raise ValueError("len(pos) must be SpaceDim")

    LoEnd,HiEnd=LD.PhysDom()
    DomBox=Box(LoEnd+HiEnd+('Box',))


    SliceBox=Box(tuple(pos)+tuple(pos)+('Box',))
    if DomBox.ContainsBox(SliceBox):
    # note that only the proc with procID=0 will have data. All other processes
    # have nans.
        try:
            return LD.GatherBox(SliceBox,components=[i for i in range(LD.nComp)])
        except:
            return LD.GatherFluxBox(SliceBox)
    else:
        return None


def PlotSlice(LD, dir=0, pos=0, comp=0,filename=None, title=''):
    """
    Plot a slice of the comp FAB or FluxBox (depending on the input) sliced over a plane with normal
    in direction dir and at logical position given by pos. If the plane is outside the
    enveloping box, returns None. Note that when the plane is within the enveloping box,
    only the root process will have actual data. If filename is None, then the plot is
    sent to the screen. Otherwise, a file is created. The type of file is determined by the extension of
    filename. See the documentation of matplotlit.pyplot.savefig for more details.
    """
    try:
        import matplotlib.pyplot as plt

    except:
        print("You do not seem to have matplotlib installed. Run pip(3) install matplotlib --user ")
        
    X=GatherSlice(LD, dir=dir, pos=pos)
    if LD.myMpiId==0:
        try:
            v=X.Fluxes[comp][...,0]
        except:
            try:
                v=X[...,comp]
            except: #the slice is empty
                v=np.empty((1,1))
        finally:
            
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            if LD.SpaceDim == 3:
                
                mesh=ax1.imshow(np.squeeze(v.transpose()), interpolation='bilinear')
                plt.colorbar(mesh)
            else:
                ax1.plot(np.squeeze(v.transpose()))
            
            plt.title(title)
            if filename is None:
                plt.show()
            else:
                plt.savefig(filename)
                plt.close()


def Plotline(LD, dir=0, pos=(0,0,0), comp=0,filename=None, title=''):
    """
    Plot a line of the comp FAB or FluxBox (depending on the input) sliced over a line
    in direction dir passing by the logical  pos in the normal plane. Note that only the values of pos
    not in direction dir are used. If the line is outside the
    enveloping box, returns None. Note that when the line is within the enveloping box,
    only the root process will have actual data.   The type of file is determined by the extension of
    filename. See the documentation of matplotlit.pyplot.savefig for more details.
    """
    
    try:
        import matplotlib.pyplot as plt
    except:
        print("You do not seem to have matplotlib installed. Run pip(3) install matplotlib --user ")

    
    X=GatherLine(LD, dir=dir, pos=pos)
    
    if LD.myMpiId==0:
        try:
            v=X.Fluxes[comp][...,0]
        except:
            try:
                v=X[...,comp]
                
            except: #the slice is empty
                v=np.empty((1,1))
        finally:
            
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
           
            ax1.plot(np.squeeze(v))

            plt.title(title)
            if filename is None:
                plt.show()
            else:
                plt.savefig(filename)
                plt.close()

    





