import numpy as np
import PyGlue as pg
#import matplotlib.pyplot as plt
#from numba import jit
import time
from ElliKit import Ellikit
from scipy.sparse import csr_matrix
from scipy.sparse import linalg as LA
import matplotlib.pyplot as plt
from SOMAR import SOMAR
M = []
X = []
bcdesc = []
Box=[]
c=0
#def GenerateMatrix(nx, ny, nz, dx, dy, dz):

def SomarToPyBox(box):
    LowCorner = [0,0,0]
    Size = [1,1,1]
    Centering = [0, 0, 0]
    SpaceDim = len(box) // 2

    for d in range(SpaceDim):
        LowCorner[d]=box[d]
        Size[d]=box[d+SpaceDim]-box[d]+1

    return (tuple(LowCorner),tuple(Size),tuple(Centering))

def SomarToPyBcDesc(bcDesc):
    R = []
    bcDescriptor = pg.ValarrayToNumpy(bcDesc)

    for bc in bcDescriptor:

        if bc==0:
            R.append(None)
        elif bc==1:
            R.append('Dirichlet')
        elif bc==2:
            R.append('Neumann')
        else:
            print(bc)
            raise ValueError("ElliKitWrapper:SomarToPyBcDesc has been called with an unrecognized bcDescriptor. \
                              This is most likely due to a box sharing a side with two or more boxes.  ")

    return tuple(R)

def compare(v1, v2):

    C1 = pg.FABToNumpy(v1)
    C2 = pg.FABToNumpy(v2)
    size = C1.shape
    print("size", size)
    index = np.arange(size[0] * size[1] * size[2])
    I = np.unravel_index(index, size[0:-1], order='F')[0]
    J = np.unravel_index(index, size[0:-1], order='F')[1]
    K=np.unravel_index(index, size, order='F')[2]
    for (i, j, k, c1, c2) in zip(I,J,K, C1.ravel('F'), C2.ravel('F')):
        if (abs(c1 - c2) > 1.e-8):
            print(i, j, k, c1, c2)

def output(v1):
    return
    C1 = pg.FABToNumpy(v1)
    size = C1.shape
    print("size", size)
    index = np.arange(size[0] * size[1] * size[2])
    I = np.unravel_index(index, size[0:-1], order='F')[0]
    J = np.unravel_index(index, size[0:-1], order='F')[1]
    K=np.unravel_index(index, size, order='F')[2]
    for (i, j, k, c1) in zip(I, J, K, C1.ravel('F')):
        print(i-1, j-1, k-1, c1)

def plotq2(v1, v2):
    C1 = pg.FABToNumpy(v1)
    C2 = pg.FABToNumpy(v2)
    print(C1.shape)
    print(C2.shape)
    c1=C1[1,:,1]
    s=c1.shape

    x1=np.arange(0,s[0])
    print(x1)
    plt.plot(x1,C1[1,:, 1])
    c2=C2[1,:,1]
    s=c2.shape
    x2=np.arange(1,s[0]+1)
    if(not c2.shape  == c1.shape):
        plt.plot(x2,C2[0,:,0])
    else:
        plt.plot(x1,C2[1,:,1])
    plt.show()

def plotq1(v1):
    C1 = pg.FABToNumpy(v1)

    c1=np.zeros((C1.shape[0],C1.shape[1]), order='F')
    for i in range(16):
        c1[:,:]=C1[:,:,0].reshape((C1.shape[0], C1.shape[1]))
        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        mesh=ax1.imshow(c1.transpose(), interpolation='bilinear')

        plt.title(i)
        plt.colorbar(mesh)
        plt.show()

    return
    s=c1.shape

    x1=np.arange(0,s[0])

    plt.plot(x1,C1[1,:, 1])

    plt.show()

def printBox(box):
    print(" *** box ***", box)
    Box=SomarToPyBox(box)
def printBC(BC):
    print(" *** BC ***", SomarToPyBcDesc(BC))


def GenerateMatrix(dx, box, refRatio, domain,BC=None, Jgup=None):



    PhiBox = SomarToPyBox(box)

    DX = dx[:-1]
    RR=refRatio[:-1]
    Domain=SomarToPyBox(domain[:-1])
    if BC is None:
        bc=None
    else:
        bc=SomarToPyBcDesc(BC)

        M.append(Ellikit(PhiBox, DX, BC=bc, Domain=Domain, RefRatio=RR))

    return

def GetMatrixnnz(OperatorName):

    Op = OperatorName[0]
    return M[0].CSRData(Op=Op)['nnz']

def GetMatrixRows(OperatorName):
    Op = OperatorName[0]

    return M[0].CSRData(Op=Op)['shape'][0]

def GetMatrixColumns(OperatorName):
    Op=OperatorName[0]
    return M[0].CSRData(Op=Op)['shape'][1]

def GetMatrixData(OperatorName, ptr, col, data):
    Op=OperatorName[0]
    Ptr = pg.ValarrayToNumpy(ptr)
    Col = pg.ValarrayToNumpy(col)
    Data = pg.ValarrayToNumpy(data)
    s=M[0].CSRData(Op=Op)['indptr'].shape[0]
    Ptr[0:s] = M[0].CSRData(Op=Op)['indptr']
    Col[:] = M[0].CSRData(Op=Op)['indices']
    Data[:] = M[0].CSRData(Op=Op)['data']

    return M[0].CSRData(Op=Op)['shape'][0]

def GetMatrixDiagonal(OperatorName, data):

    Op = OperatorName[0]
    data = pg.ValarrayToNumpy(data)
    data[:] = M[0].GetDiagonals(Op=Op)

def ClearMatrix():
    del M[0]

def PySolution(x):

    R = pg.ValarrayToNumpy(x)
    b = np.ones(M[0]._Laplacian[0].shape[0], dtype='double')
    start = time.time()

    R[:] = LA.cg(M[0]._Laplacian[0], b, tol=1.e-5, callback=counter, maxiter=10000)[0]
    print("time for Python solution is ", time.time()-start, "s")
    print("Error of Python solution is ", np.sqrt(np.linalg.norm(M[0]._Laplacian[0].dot(R) - b)), " after ", c, "iterations")

    return

def counter(x):
    global c
    c+=1
    return c

def LowestMode(col, rowOffset, Data, eig):
    print(" in LowestMode")
    columns=pg.ValarrayToNumpy(col)
    row_offsets=pg.ValarrayToNumpy(rowOffset)
    values=pg.ValarrayToNumpy(Data)
    EIG=pg.ValarrayToNumpy(eig)
    nrows=row_offsets.shape[0]-1
    print(" nrows " ,nrows)
    M=csr_matrix((values,columns, row_offsets), shape=(nrows, nrows))

    eigenvalues, eigenvectors = LA.eigsh(M, k=2, which='SM')
    EIG[:]=eigenvectors[:,0]
    print(" eigen values ", eigenvalues)

    return














