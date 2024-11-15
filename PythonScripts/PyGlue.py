from mpi4py import MPI
import numpy as np
comm = MPI.COMM_WORLD
MpiNumProcs = comm.Get_size()
myMpiId = comm.Get_rank()
valarrayTypes = {'int': np.int32,'float': np.float32,'double': np.float64, 'size_t' : np.uint64}


def make_from_memView(memView, shape, dtype=np.float64, order='C', own_data=True):
    """ creates an numpy array of shape shape and type dtype
        which points to an address given by memView. The order can be C (last index runs first)
        or F (first index runs first)"""

    arr = np.ndarray(tuple(shape[:]), dtype, memView, order=order)
    return arr

def FABToNumpy(Args):
    """unpacks the FAB specified by Args into a numpy: Args[0]-> tuple of lower index of box
                                                       Args[1]-> tuple of size
                                                       Args[2]-> tuple of centering
                                                       Args[3]-> number of components
                                                       Args[4]-> address of buffer """
    size=Args[1]
    nComp=Args[3]
    view=Args[4]
    arr = make_from_memView(view,np.append(size,np.array(nComp)),order='F')
    return arr

def ValarrayToNumpy(v):
    """unpacks valarray (or array) specified by v"""
    size=v[0]
    TypeOfObject=v[1]
    view=v[2]

    arr=make_from_memView(view,[size],order='F', dtype=valarrayTypes[TypeOfObject])
    return arr

def BoxOfFAB(Args):
    """unpacks the box of a FAB specified by Args. On output returs a tuple of tuples containing
        lower corner, higher corner, centering"""
    box = Args[0]
    size = Args[1]
    hi=[]
    for i in range(len(box)):
        hi.append(box[i] + size[i]-1)


    box += tuple(hi)

    centering = Args[2]
    box += centering
    return box


def FABToCoordinates(Args):
    Coord = FABToNumpy(Args)
    if len(Coord.shape)==4: # 3D case
        return Coord[..., 0], Coord[..., 1], Coord[..., 2]
    else:
        return Coord[..., 0], Coord[..., 1]



def CoordToIndex():
    pass







def PhysDom(Arg):
    """ Determines the extent of the domain from a leveldata. Note that this may not the total extent of the
    actual domain. It is just the smallest box that contains all the boxes in the LD collection """



    lo = [x for x in Arg[1][0][0]]
    hi = [x for x in Arg[1][0][1]]

    for Box in Arg[1]:
        for d in range(len(lo)):
            lo[d] = min(lo[d], Box[0][d])
            hi[d] = max(hi[d], Box[1][d])

    return (tuple(lo), tuple(hi))
