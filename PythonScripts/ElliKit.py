import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, dia_matrix, coo_matrix, identity, dok_matrix
from scipy.sparse import linalg as LA
import time
#import sys
#import warnings
#if not sys.warnoptions:
 #   warnings.simplefilter("ignore")
    # removes warning about not using lil_matrix format
#global dir3D
dir3D = ((1, 0, 0), (0, 1, 0), (0, 0, 1))


class Ellikit:
    def __init__(self,box, DX, Jgup=None, Jinv=None, BC=None, Domain=None, RefRatio=None):
        """Produces the suite of vector calculus operators applicable to fields defined in box.
        For each direction, two boundary conditions are required. The BC are None, Dirichlet or Neumann.
        If None, the operator assumes that a ghost layer in the corresponding direction  and sidedeness
        is present on the fields it acts upon. Otherwise, the ghost layer is virtually generated assuming an extrapolation from
        the interior via mirror symmetry (Neumann) or antisymmetry (Dirichlet). In this latter cases, the field must not have a ghost layer.
        In all cases, the operators return a field defined on the interior points (see examples at the bottom) ."""
        self.time=dict()

        self.GridLowerCorner = box[0]
        self.GridSize = box[1]
        self.GridCentering = box[2]
        if Domain is None: # assume Domain=box
            self.DomainLowerCorner = self.GridLowerCorner
            self.DomainSize = self.GridSize

        else:
            self.DomainLowerCorner = Domain[0]
            self.DomainSize = Domain[1]

        self.DX = DX
        self.NDim = len([i for i in self.GridSize if i > 1])
        if BC is None:
            self.BC=tuple([None for i in range(2*self.NDim)])
        else:
            self.BC=BC

        self.RefRatio=RefRatio
        # private members have leading underscores (so to speak, Python has no real concept of "private")
        self._Grad = []  # Grad
        self._setGrad()
        self._Div=[]
        self._setDiv() # Div is defined on the flux box that wraps the interior valid points
        self._Laplacian = self._laplacian()
        self._Outer, self._Inner = self._InnerOuterSplit()
        self._LaplacianNoGhosts=self._deGhost()
        self._LaplacianOnGhosts=self._onGhost()
        if RefRatio is not None:
            self._Restriction = self._restriction()
            #print("Restriction operator generated with RefRatio ", RefRatio)
            self._Prolongation = self._prolongation()
        self._MainDiagonals = self._diagonals(self._LaplacianNoGhosts[0:-1])
        #self._D2G = self._DomainToGrid(self.DomainLowerCorner, self.DomainSize, self.GridLowerCorner, self.GridSize)
        #calculates the reduced grid. This is the grid reduced along boundaries where the Boundary condition is None
        GridLC = [0,0,0]
        GridS = [1, 1, 1]
        DomainLC = [0, 0, 0]
        DomainS = [1, 1, 1]

        for d, bc in enumerate(self.BC):

            GridLC[d // 2] = self.GridLowerCorner[d // 2]
            GridS[d // 2] = self.GridSize[d // 2]
            DomainLC[d // 2] = self.DomainLowerCorner[d // 2]
            DomainS[d //2] = self.DomainSize[d // 2]

        for d, bc in enumerate(self.BC):
            if bc is None and d % 2 == 0:
                 # lower boundary
                GridLC[d // 2] += 1
                GridS[d // 2] -= 1
            if bc is None and d % 2 == 1:
                #upper boundary
                GridS[d // 2] -= 1

        self._G2D = (self._GridToDomain(DomainLC, DomainS, GridLC, GridS),GridS)


    # public operators






    def CSRData(self, Op=None):
        """ Returns the components and other details of  the Operator Op in CSR format. Returns a dictionary """
        if Op=='Laplacian':
            return self._returnCSRData(self._Laplacian[0])

        elif Op == 'LaplacianDMinusOne':
            return self._returnCSRData(self._Laplacian[1])

        elif Op == 'Restriction':
            return self._returnCSRData(self._Restriction[0])
        elif Op == 'Prolongation':
            return self._returnCSRData(self._Prolongation[0])
        elif Op== 'LaplacianNG':
            return self._returnCSRData(self._LaplacianNoGhosts[0])
        elif Op== 'LaplacianNGDMinusOne':
            return self._returnCSRData(self._LaplacianNoGhosts[1])
        elif Op== 'LaplacianOG':
            return self._returnCSRData(self._LaplacianOnGhosts[0])
        elif Op== 'LaplacianOGDMinusOne':
            return self._returnCSRData(self._LaplacianOnGhosts[1])
        elif Op == 'Divergence x':
            return self._returnCSRData(self._Div[0])
        elif Op == 'Divergence y':
            if self.NDim < 2:
                raise ValueError("Ellikit: CSRData was defined with in less than 2D, but Divergence y was requested.")
            return self._returnCSRData(self._Div[1])
        elif Op == 'Divergence z':
            if self.NDim < 3:
                raise ValueError("Ellikit: CSRData was defined with in less than 3D, but Divergence z was requested.")

            return self._returnCSRData(self._Div[2])
        elif Op == 'Gradient x':

            return self._returnCSRData(self._Grad[0])
        elif Op == 'Gradient y':
            return self._returnCSRData(self._Grad[1])
        elif Op == 'Gradient z':
            if self.NDim < 3:
                raise ValueError("Ellikit: CSRData was defined with in less than 3D, but Gradient z was requested.")

            return self._returnCSRData(self._Grad[2])
        elif Op == 'DomainToGrid':
            return self._returnCSRData(self._D2G)
        elif Op == 'GridToDomain':
            return self._returnCSRData(self._G2D[0])
        else:
            raise ValueError("ElliKit: CSRData called with the wrong/not implemented yet operator")

    def GetDiagonals(self, Op=None):

        if Op == 'Laplacian' or Op=='LaplacianNG':
            return self._MainDiagonals[0]
        elif Op=='LaplacianDMinusOne' or Op=='LaplacianNGDMinusOne':
            return self._MainDiagonals[1]
        else:
            raise ValueError("ElliKit: Diagonals called with the wrong/not implemented yet operator")


    def Dot(self,v, Op=None):
        """Applies the operator Op to the field v """
        if Op == 'Divergence':
            r = np.zeros(self._Div[0].shape[0])
            for d,q in zip(self._Div,v):
                r += d.dot(q.ravel('F'))
            r.reshape(self._Laplacian[2], order='F')
            return r
        elif Op=='Laplacian':
            for i, o in zip(v.shape, self.GridSize):
                if i != o:
                    raise ValueError("Ellikit Dot: the shape of the input does not match the shape the operator was defined")

            start=time.time()
            r = self._Laplacian[0].dot(v.ravel('F')).reshape(self._Laplacian[2], order='F')
            self._addTime('LaplacianDot', time.time() - start)
            return r
        elif Op=='LaplacianDMinusOne':
            for i, o in zip(v.shape, self.GridSize):
                if i != o:
                    raise ValueError("Ellikit Dot: the shape of the input does not match the shape the operator was defined")

            start=time.time()
            r = self._Laplacian[1].dot(v.ravel('F')).reshape(self._Laplacian[2], order='F')
            self._addTime('LaplacianDMinusOneDot', time.time()-start)
            return r

        elif Op == 'LaplacianNG':
            for i, o in zip(v.shape, self._LaplacianNoGhosts[2]):
                if i != o:
                    raise ValueError("Ellikit Dot: the shape of the input does not match the shape the operator was defined")

            start=time.time()
            r = self._LaplacianNoGhosts[0].dot(v.ravel('F')).reshape(self._LaplacianNoGhosts[2], order='F')
            self._addTime('LaplacianNGDot', time.time() - start)
            return r
        elif Op == 'LaplacianNGDMinusOne':
            for i, o in zip(v.shape, self._LaplacianNoGhosts[2]):
                if i != o:
                    raise ValueError("Ellikit Dot: the shape of the input does not match the shape the operator was defined")

            start=time.time()
            r = self._LaplacianNoGhosts[1].dot(v.ravel('F')).reshape(self._LaplacianNoGhosts[2], order='F')
            self._addTime('LaplacianNGDMinusOneDot', time.time() - start)
            return r
        elif Op == 'LaplacianOG':
            for i, o in zip(v.shape, self.GridSize):
                if i != o:
                    raise ValueError("Ellikit Dot: the shape of the input does not match the shape the operator was defined")

            start=time.time()
            r = self._LaplacianOnGhosts[0].dot(v.ravel('F')).reshape(self._LaplacianOnGhosts[2], order='F')
            self._addTime('LaplacianOGDot', time.time() - start)
            return r
        elif Op == 'LaplacianOGDMinusOne':
            for i, o in zip(v.shape, self.GridSize):
                if i != o:
                    raise ValueError("Ellikit Dot: the shape of the input does not match the shape the operator was defined")

            start=time.time()
            r = self._LaplacianOnGhosts[1].dot(v.ravel('F')).reshape(self._LaplacianOnGhosts[2], order='F')
            self._addTime('LaplacianOGDot', time.time() - start)
            return r
        elif Op == 'Inner':
            for i, o in zip(v.shape, self.GridSize):
                if i != o:
                    raise ValueError("Ellikit Dot: the shape of the input does not match the shape the operator was defined")

            start=time.time()
            r = self._Inner[0].dot(v.ravel('F')).reshape(self._Inner[1], order='F')
            self._addTime('InnerDot', time.time() - start)
            return r

        elif Op == 'Outer':
            for i, o in zip(v.shape, self.GridSize):
                if i != o:
                    raise ValueError("Ellikit Dot: the shape of the input does not match the shape the operator was defined")

            start=time.time()
            r = self._Outer[0].dot(v.ravel('F')).reshape(self._Outer[1], order='F')
            self._addTime('OuterDot', time.time() - start)
            return r
        elif Op=='Restriction':
            for i, o in zip(v.shape, self.GridSize):
                if i != o:
                    raise ValueError("Ellikit Dot: the shape of the input does not match the shape the operator was defined")

            start=time.time()
            GS=self.GridSize
            RR=self.RefRatio
            shape=(GS[0]//RR[0], GS[1]//RR[1], GS[2]//RR[2])
            r= self._Restriction[0].dot(v.ravel('F')).reshape(shape, order='F')
            self._addTime('Restriction', time.time()-start)
            return r
        else:
            raise ValueError("ElliKit: Dot called with the wrong/not implemented yet operator")



    # private operators
    def _DomainToGrid(self, DomainLC, DomainS, GridLC, GridS):
        """ Returns the Matrix representation of the operator that takes a field defined over the whole domain
        and returns that is supported by self.Grid """
        if len(GridLC) == 1:
            return self._DomainToGrid((DomainLC[0],0,0), (DomainS[0],1,1), (GridLC[0],0,0), (GridS[0], 1, 1))
        if len(GridLC) == 2:
            return self._DomainToGrid((DomainLC[0], DomainLC[1], 0), (DomainS[0], DomainS[1], 1),
            (GridLC[0],GridLC[1],0), (GridS[0], GridS[1], 1))
        start=time.time()
        Nc = 1
        Nr = 1

        for d,g in zip(DomainS, GridS):
            Nc *= d
            Nr *= g

        diag = np.zeros(Nc, dtype='int')
        indicesR = np.arange(Nc).reshape(DomainS, order='F')
        indicesR.shape
        GridLowerCorner=np.zeros(3, dtype='int')
        for n,(g, d) in enumerate(zip(GridLC, DomainLC)):
            if (g - d < 0) or ((g-d)+GridS[n] > DomainS[n]) :
                raise ValueError("ElliKit: DomainToGrid has been called with the grid outside the domain")
            GridLowerCorner[n]=g - d

        ind = indicesR[GridLowerCorner[0]:GridLowerCorner[0]+GridS[0],
                       GridLowerCorner[1]:GridLowerCorner[1]+GridS[1],
                       GridLowerCorner[2]:GridLowerCorner[2]+GridS[2]].ravel('F')
        diag[ind] = 1

        R = dia_matrix((diag, np.zeros(1, dtype='int')), shape=(Nc, Nc)).tocsr()
        mask = np.concatenate(([True], R.indptr[1:] != R.indptr[:-1]))

        R = csr_matrix((R.data, R.indices, R.indptr[mask]), shape=(Nr, Nc))

        self._addTime('DomainToGrid', time.time() - start)
        return R





    def _GridToDomain(self, DomainLC, DomainS, GridLC, GridS):
        """ Returns the Matrix representation of the operator that takes a field defined over the grid
        and extends it by zero padding to the whole domain. It is calculated taking the transpose of the correspondig
        DomainToGrid operator """
        start=time.time()
        if len(GridLC) == 1:
            return self._DomainToGrid((DomainLC[0],0,0), (DomainS[0],1,1), (GridLC[0],0,0), (GridS[0], 1, 1))
        if len(GridLC) == 2:
            return self._DomainToGrid((DomainLC[0], DomainLC[1], 0), (DomainS[0], DomainS[1], 1),
            (GridLC[0],GridLC[1],0), (GridS[0], GridS[1], 1))
        start=time.time()
        Nc = 1
        Nr = 1

        for d,g in zip(DomainS, GridS):
            Nc *= g
            Nr *= d

        diag = np.zeros(Nr, dtype='int')
        indicesC = np.arange(Nr).reshape(DomainS, order='F')

        GridLowerCorner=np.zeros(3, dtype='int')
        for n,(g, d) in enumerate(zip(GridLC, DomainLC)):
            if (g - d < 0) or ((g-d)+GridS[n] > DomainS[n]) :
                raise ValueError("ElliKit: DomainToGrid has been called with the grid outside the domain")
            GridLowerCorner[n]=g - d

        ind = indicesC[GridLowerCorner[0]:GridLowerCorner[0]+GridS[0],
                       GridLowerCorner[1]:GridLowerCorner[1]+GridS[1],
                       GridLowerCorner[2]:GridLowerCorner[2]+GridS[2]].ravel('F')
        diag[ind] = 1

        R = dia_matrix((diag, np.zeros(1, dtype='int')), shape=(Nr, Nr)).tocsc()
        mask = np.concatenate(([True], R.indptr[1:] != R.indptr[:-1]))

        R = csc_matrix((R.data, R.indices, R.indptr[mask]), shape=(Nr, Nc)).tocsr()


        self._addTime('GridToDomain', time.time() - start)
        return R



    def _deGhost(self):
        """ Combines the Laplacian(s) with the Inner matrix, and eliminates the zero columns,
        returning a square matrix that operates on the interior points"""
        R=[]
        for M in self._Laplacian[0:-1]:
            if M is not None:
                NG=M.dot(self._Inner[0])
                R.append(self._stripOuterColumns(NG))
        Ns=self.GridSize
        BC=[(i is None) for i in self.BC]
        L=[int(i) for (n,i) in enumerate(BC) if n%2==0] # lower limits (0 if BC is not None, 1 if it is)
        U=[Ns[(n-1)//2]-int(i) for (n,i) in enumerate(BC) if n%2==1] # upper limits (gotta love list comprehension)
        NInner=tuple([u-l for u,l in zip(U,L)])
        R.append(NInner)
        return tuple(R)

    def _onGhost(self):
        """ Combine the Laplacians(s) with the Outer matrix and eliminates the zero entries."""
        R=[]
        for M in self._Laplacian[0:-1]:
            if M is not None:
                OG=M.dot(self._Outer[0])
                OG.eliminate_zeros()
                R.append(OG)
        Ns=self.GridSize
        BC=[(i is None) for i in self.BC]
        L=[int(i) for (n,i) in enumerate(BC) if n%2==0] # lower limits (0 if BC is not None, 1 if it is)
        U=[Ns[(n-1)//2]-int(i) for (n,i) in enumerate(BC) if n%2==1] # upper limits (gotta love list comprehension)
        NInner=tuple([u-l for u,l in zip(U,L)])
        R.append(NInner)
        return tuple(R)

    def _stripOuterColumns(self, A):
        """Helper function  that eliminates the columns corresponding to outer points, as
        defined by the BCs. It assumes that those columns have been previously set to zero
        by multiplication with the Inner matrix. """
        Ns=self.GridSize
        A.eliminate_zeros() # The columns corresponding to ghost points are zero
        A.sum_duplicates() # just to be sure
        if  not A.has_sorted_indices:
            A.sort_indices()

        # all we need to do is to rebase the column index
        BC=[(i is None) for i in self.BC]

        L=[int(i) for (n,i) in enumerate(BC) if n%2==0] # lower limits (0 if BC is not None, 1 if it is)
        U=[Ns[(n-1)//2]-int(i) for (n,i) in enumerate(BC) if n%2==1] # upper limits (gotta love list comprehension)

        NInner=tuple([u-l for u,l in zip(U,L)])

        C1=tuple([c1-l for (n,l),c1 in zip(enumerate(L),np.unravel_index(A.indices, Ns,order='F'))])
        newcol=np.ravel_multi_index(C1,NInner,order='F')

        return csr_matrix((A.data, newcol, A.indptr), (A.shape[0], A.shape[0]))

    def _returnCSRData(self, Matrix):
        """ used to create a dictionary with the data in a CSRMatrix"""
        return {'data':Matrix.data, 'indices':Matrix.indices,
                'indptr': Matrix.indptr, 'nnz': Matrix.nnz,
                'shape': Matrix.shape, 'AreIndicesSorted': Matrix.has_sorted_indices}


    def _addTime(self, key, t):
        """ Accumulates t into self.time['key'] """
        if key in self.time:
            self.time[key] += t
        else:
            self.time[key] = t

    def _Split(self, dir, Side, Ns, *args):
        """ Generates a pair of matrices which separates the boundary from the interior.
        The sum is the identity matrix. If args[0]=None, it simply returns the identity matrix"""
        # sanity checks

        if dir not in range(len(Ns)):
            raise ValueError(" Ellikit.Split: dir must be < len(Ns)")
        if Side not in (0, 1):
            raise ValueError("Ellikit.Split: Side must be 0 or 1")
        if len(Ns) == 1:
            return self._Split(dir, Side, (Ns[0], 1, 1), *args)
        if len(Ns) == 2:
            return self._Split(dir, Side, (Ns[0], Ns[1], 1), *args)
        start = time.time()
        N = 1
        for i in Ns:
            N *= i

        if len(args) > 0:
            if args[0] is None:
                return csr_matrix((N,N), dtype='int'), identity(N, dtype='int', format='csr')
            else:
                raise ValueError("Ellikit.Split: called with an extra argument which was not None")



        indicesR = np.arange(N).reshape(Ns, order='F')
        if dir == 0:
            ind = indicesR[-Side,:,:].ravel('F')
        elif dir == 1:
            ind = indicesR[:, -Side,:].ravel('F')
        elif dir ==2:
            ind = indicesR[:,:, -Side].ravel('F')
        diag = np.ones(N, dtype='int')
        diag[ind] = 0

        Inner = dia_matrix((diag, np.zeros(1, dtype='int')), shape=(N,N))

        Outer = identity(N, dtype='int', format='dia')
        Outer -=  Inner
        self._addTime('Split', time.time()-start)
        return (Outer.tocsr(),Inner.tocsr())

    def _AddGhostLayer(self, dir, Side, Ns, TYPE='None'):
        """Returns a CSR matrix that adds a ghost layer. The input field has dimensions Ns, and
        the output field has the dir dimension increased by one if the TYPE of BC is other than None.
        Side determines if the ghost layer
        is added the beginning or end. If BC is specified, the ghost layer is set to minus its mirror image
        if BC = 'Dirichlet', or to its mirror image if BC = 'Neumann'. """
        # sanity checks
        if TYPE not in (None, 'Dirichlet', 'Neumann'):
            raise ValueError(
                " Ellikit.AddGhostLayer: BC must be one of these: None, 'Dirichlet' or 'Neumann'")
        if dir not in range(len(Ns)):
            raise ValueError(" Ellikit.AddGhostLayer: dir must be < len(Ns)")
        if Side not in (0, 1):
            raise ValueError("Ellikit.AddGhostLayer: Side must be 0 or 1")
        if len(Ns) == 1:
            return self._AddGhostLayer(dir, Side, (Ns[0], 1, 1), TYPE=TYPE)
        if len(Ns) == 2:
            return self._AddGhostLayer(dir, Side, (Ns[0], Ns[1], 1), TYPE=TYPE)

        start = time.time()
        N, n = 1, 1
        for i, s in zip(Ns, dir3D[dir]):
            if i>1:
                N *= i
                n *= (i + s)
        if TYPE is None:
            return identity(N, dtype='int', format='csr')
        else:
            s = dir3D[dir]
            Side = 1-Side  # we will undo that later
            L = (Side * s[0], Side * s[1], Side * s[2])
            H = (Ns[0]+Side*s[0], Ns[1]+Side*s[1], Ns[2]+Side*s[2])
            indicesC = np.arange(N).reshape(Ns, order='F').ravel('F') # column indices
            indicesR = np.arange(n).reshape(
                (Ns[0] + s[0], Ns[1] + s[1], Ns[2] + s[2]), order='F') # row indices of the output
            ind = indicesR[L[0]:H[0], L[1]:H[1], L[2]:H[2]].ravel('F')  # rows corresponding to the "interior"
            # i.e., the identity matrix is given by (ind, indicesC)

            Side = 1-Side

            # below, ind2 is the row that provides the values that are extrapolated to the row ind1
            if dir == 0:
                if Side == 0:
                    ind1 = indicesR[0, :, :].ravel('F') # row that needs to be extrapolated
                    ind2 = indicesR[1, :, :].ravel('F') # row that provides the material for extrapolation.

                elif Side == 1:
                    ind1 = indicesR[-1, :, :].ravel('F')
                    ind2 = indicesR[-2, :, :].ravel('F')

            elif dir == 1:
                if Side == 0:
                    ind1 = indicesR[:, 0, :].ravel('F')
                    ind2 = indicesR[:, 1, :].ravel('F')

                elif Side == 1:
                    ind1 = indicesR[:, -1, :].ravel('F')
                    ind2 = indicesR[:, -2, :].ravel('F')
            elif dir == 2:
                if Side == 0:
                    ind1 = indicesR[:, :, 0].ravel('F')
                    ind2 = indicesR[:, :, 1].ravel('F')
                elif Side == 1:
                    ind1 = indicesR[:, :, -1].ravel('F')
                    ind2 = indicesR[:, :, -2].ravel('F')

            nn = ind1.size

            if TYPE == 'Dirichlet':
                data1=-np.ones(nn, dtype='int')
            elif TYPE == 'Neumann':
                data1 = np.ones(nn, dtype='int')
            ind1C = np.empty_like(data1, dtype='int')
            # now we fill ind1C with the corresponding column value
            # first we create a lookup table, that gives the row
            # that corresponds to a given position along ind1.
            # So for a given value, we find its position in ind
        ind_dict = {value: pos for pos, value in enumerate(ind)}
            # next we go over all values  in ind2, and we find the
            # corresponding position in ind1 via the lookup table.
            # in indicesC, at that position, we find the column.
        for (i,_),j in zip(enumerate(ind1), ind2):
            pos = ind_dict[j]
            ind1C[i] = indicesC[pos]
        # Finally, we add the values to data
        data = np.concatenate((np.ones_like(ind, dtype='double'), data1))
        ind = np.concatenate((ind, ind1)) # add the ind1 of the top or bottom row
        indicesC = np.concatenate((indicesC, ind1C)) # add the columns.
        M = csr_matrix(coo_matrix(
            (data, (ind, indicesC)), shape=(n, N)))
        self._addTime('AddGhostLayer', time.time()-start)

        return M

    def _diagonals(self, Matrices, Ns=None, BC=None):
        """ Extracts the main diagonal from the Laplacian with BC=None. This requires some calculations since
        it is not a square matrix. For a Laplacian where all the BC are not None, it is a square matrix and it is easy"""
        if Ns is None:
            Ns=[i for i in self.GridSize if i>1]
        if BC is None:
            BC=self.BC
        if len(Ns)==1:
            return self._diagonals(Matrices, Ns=(*Ns,1,1),BC=(*BC, 'X', 'X', 'X', 'X'))
        elif len(Ns)==2:
            return self._diagonals(Matrices, Ns=(*Ns,1), BC=(*BC, 'X', 'X'))
        if 2*len(Ns) != len(BC):
            raise ValueError(" There  must be two BCs for each direction in Ns")

        start=time.time()
        R = []
        BC=[(i is None) for i in BC]

        L=[int(i) for (n,i) in enumerate(BC) if n%2==0] # lower limits (0 if BC is not None, 1 if it is)
        U=[Ns[(n-1)//2]-int(i) for (n,i) in enumerate(BC) if n%2==1] # upper limits (gotta love list comprehension)
        NInner=tuple([u-l for u,l in zip(U,L)])
        for M in Matrices:


            if M.shape[0] == M.shape[1]:  # square matrix
                R.append(M.diagonal())
            else:
                A=M.dot(self._Inner[0]) # the Laplacian localized on the inner points
                A.eliminate_zeros() # The columns corresponding to ghost points are zero
                A.sum_duplicates() # just to be sure
                if  not A.has_sorted_indices:
                    A.sort_indices()

                # all we need to do is to rebase the column index
                col=A.indices
                C1=list(np.unravel_index(col, Ns,order='F'))

                for n,l in enumerate(L):
                    C1[n]-=l

                C1=tuple(C1)


                newcol=np.ravel_multi_index(C1,NInner,order='F')
                R.append(csr_matrix((A.data,newcol, A.indptr), (M.shape[0], M.shape[0])).diagonal())


        self._addTime("LaplaceDiagonals", time.time() - start)
        return tuple(R)





    def _filter_list(self, full_list, excludes):
        """returns a lazy iterable that returns the elements in full_list that are not in S """
        S = set(excludes)
        return (x for x in full_list if x not in S)

    def _setDiv(self):
        # sanity checks
        for dir in range(self.NDim):

            # we only define it for a fluxBox created over the interior (valid) points
            Ns=list(self.GridSize)
            for d, i in enumerate(Ns):
                if i > 1:
                    Ns[d] += -(self.BC[2 * d] is None) - (self.BC[2 * d + 1] is None) + dir3D[dir][d]

            self._Div.append(self._grad(dir, Ns, self.DX[dir]))
        return

    def _grad(self, dir, Ns, delta):
        """ Returns gradient in direction dir. Ns is a tuple of dimensions of the input array.
        The output array has the dimension in the direction dir reduced by one. Assumes Fortran
        ordering of arrays.The spacing delta is along the dir direction """
        # sanity checks
        if dir not in range(len(Ns)):
            raise ValueError(" Ellikit.AddGhostLayer: dir must be 0,1 or 2")
        if len(Ns) == 1:
            return self._grad(dir, (Ns[0], 1, 1), delta)
        if len(Ns) == 2:
            return self._grad(dir, (Ns[0], Ns[1], 1), delta)

        start = time.time()
        N, n = 1, 1
        for i, s in zip(Ns, dir3D[dir]):
            N *= i
            n *= (i - s)

        offsets = [1, Ns[0], Ns[0]*Ns[1]]
        diags = np.empty((2, N))
        diags[0, :] = np.full(N, -1./delta)
        diags[1, :] = np.full(N, 1./delta)
        s = dir3D[dir]
        M = csr_matrix(dia_matrix((diags, [0, offsets[dir]]), shape=(N, N)))
        indices = np.arange(N).reshape(Ns, order='F')[0:(Ns[0] - s[0]), 0:(Ns[1] - s[1]), 0:(Ns[2] - s[2])].ravel('F')
        M = M[indices, :]

        self._addTime('Grad', time.time()-start)

        return M

    def _InnerOuterSplit(self, Ns=None, BC=None):
        """ Calculates two matrices that together partition a field into inner and outer
        (i.e. only ghost along the internal boundaries) parts.
        The sum of outer and inner is the indentity. BC is a list of bools: True if side is internal,
        False otherwise. If not specified, it is deduced from the BC's of the class, where None denotes
        and internal boundary, and Dirichlet or Neumann a physical boundary. """
        if Ns is None:
            Ns=[i for i in self.GridSize if i>1]
        if BC is None:
            BC=[(i is None) for i in self.BC]
        if len(BC) < len([i for i in Ns if i > 1]) * 2:
            raise ValueError(
                        "Ellikit InnerOuterSplit:  BC must have twice the size of Ns")
        if len(Ns) == 1:
            return self._InnerOuterSplit((Ns[0], 1, 1),  (*BC, False, False, False, False))
        if len(Ns) == 2:
            return self._InnerOuterSplit((Ns[0], Ns[1], 1),  (*BC, False, False))

        start = time.time()

        L=[int(i) for (n,i) in enumerate(BC) if n%2==0] # lower limits (0 if BC is not None, 1 if it is)
        U=[Ns[(n-1)//2]-int(i) for (n,i) in enumerate(BC) if n%2==1] # upper limits (gotta love list comprehension)

        N=1
            # determine total size of matrix
        for i in Ns:
            N *= i

        data=np.ones(N, dtype='int')
        InternalIndices=np.arange(N).reshape(Ns, order='F')[L[0]:U[0], L[1]:U[1], L[2]:U[2]].ravel('F')
        data[InternalIndices]=0
        dataOuter=data
        Outer=csr_matrix(dia_matrix((dataOuter, np.zeros(1, dtype='int')), shape=(N,N)))
        Inner = identity(N, dtype='int', format='csr')-Outer #

        self._addTime("InnerOuterSplit", time.time()-start)
        return ((Outer, Ns),  (Inner, Ns))


    def _setGrad(self):
        Ns = []
        delta = []
        BC=[]
        Ns[:] = self.GridSize
        delta[:] = self.DX
        BC[:] = self.BC
        for d in range(3 - len(Ns)):
            Ns.append(1)

        for B in BC:
            if B not in ('Dirichlet', 'Neumann', None):
                raise ValueError(
                    "Ellikit Gradient:  BC must be one of these: None, Dirichlet or Neumann")

        start = time.time()
        ActualDimension=len([i for i in Ns if i > 1])
        for dir in range(ActualDimension):
            pos = 2*dir
            Nr = list(Ns)  # Nr keeps tab of the # of rows in the matrices

            BCL = self._AddGhostLayer(dir, 0, Nr, TYPE=BC[pos])
            Nr[dir] += BC[pos] is not None  # add one if a ghost layer was added
            BCR = self._AddGhostLayer(dir, 1, Nr, TYPE=BC[pos + 1]).dot(BCL)
            del BCL
            Nr[dir] += BC[pos+1] is not None # add one if a ghost layer was added
            #
            GL = self._grad(dir, Nr, delta[dir]).dot(BCR)
            del BCR
            self._Grad.append(GL)  # store grad including BCs
            Nr[dir]+= -1
            BCtemp = list(BC)
            BCtemp[2 * dir] = 'Something' # Anything but None will do.
            BCtemp[2 * dir + 1] = 'Something'
            indForGrad = self._stripInactiveRows(dir, Nr, BCtemp)  # we calculate here, to be used later to eliminate the
            self._Grad[dir]=self._Grad[dir][indForGrad,:]                                                        #inactive rows of the gradient.



        self._addTime('setGrad', time.time() - start)

        return


    def _laplacian(self, Ns=None, delta=None, BC=None):
        """ calculates the CSR of a laplacian. Ns is the dimension of the array, delta the grid spacing. BC
        tells what boundary conditions to use along which border (lowside, highside).
        If BC=None, then the array is assumed to have ghost layers. Otherwise, 2*Dim BCs are expected, and the
        operator now will work on fields without ghost layers.  The Laplacian is constructed adding Dim blocks,
        each of type GRxGLxBCRxBCL, with BCL(R) the matrices that adds the low (high) side Ghost layer, GL
        the contravariant gradient and GR the covariant gradient normalized by J.
        """
        if Ns is None:
            Ns=[i for i in self.GridSize if i>1]
        if delta is None:
            delta=self.DX
        if BC is None:
            BC=self.BC

        if len(Ns) == 1:
            return self._laplacian((Ns[0], 1, 1), delta, BC=BC)
        if len(Ns) == 2:
            return self._laplacian((Ns[0], Ns[1], 1), delta, BC=BC)
        # if all((i == None) for i in BC):
        #     return self._laplacianNoBC(Ns, delta) # laplacianNoBC does not define the gradient properly!!

        # sanity check
        for B in BC:
            if B not in ('Dirichlet', 'Neumann', None):
                raise ValueError(
                    "Ellikit Laplacian:  BC must be one of these: None, Dirichlet or Neumann")

        start = time.time()

        Nr,Nc = 1,1
        # determine total size of Laplacian
        for dir, i in enumerate(Ns):
            if i>1:
                pos = 2 * dir
                Nr *= i -2+ (BC[pos] is not None) + (BC[pos + 1] is not None)  #
                Nc*=i
        n = []

        Laplacian = csr_matrix((Nr, Nc))  # empty matrix to accumulate Laplacian
        ActualDimension=len([i for i in Ns if i > 1])
        LaplacianDMinusOne = None

        for dir in range(ActualDimension):
            pos = 2*dir
            Nr = list(Ns)  # Nr keeps tab of the # of rows in the matrices
            Nr[dir] += BC[pos] is not None  # add one if a ghost layer was added to the left
            Nr[dir] += BC[pos+1] is not None # add one if a ghost layer was added to the right
            #
            Nr[dir]+= -1 # gradient and divergence each take one out.
            Nr[dir]+= -1
            Laplacian += self._Div[dir].dot(self._Grad[dir])  # accumulates in Laplacian
            n.append(Nr[dir])
            if (dir == ActualDimension - 2):
                LaplacianDMinusOne = csr_matrix(Laplacian)

        self._addTime('LaplacianConst', time.time() - start)
        return (Laplacian,LaplacianDMinusOne, tuple(n))




    def _stripInactiveRows(self,dir, N, BC=[None]*6):
        """ Eliminates rows corresponding to ghost layers in the Laplacian output"""
        Dim = len([i for i in N if i > 1])

        indices = np.arange(N[0] * N[1] * N[2]).reshape(N, order='F')
        l = [0, 0, 0]
        h = list(N)

        for d, i in enumerate(N):
            if i > 1:
                l[d] +=(BC[2*d] is  None)
                h[d] +=- (BC[2*d+1] is  None)

        if Dim == 1:  # we do not need to do anything
            return indices.ravel('F')
        elif Dim == 2:
            if dir == 0:
                return indices[:,l[1]:h[1] ,:].ravel('F')
            else:
                return indices[l[0]:h[0],:,:].ravel('F')
        else: #Dim =3
            if dir == 0:
                return indices[:,l[1]:h[1],l[2]:h[2]].ravel('F')
            elif dir == 1:
                return indices[l[0]:h[0],:,l[2]:h[2]].ravel('F')
            else:
                return indices[l[0]:h[0],l[1]:h[1],:].ravel('F')

    def _laplacianNoBC(self, Ns, delta):
        """ calculates the CSR of a laplacian. Ns is the dimension of the array (with ghosts), delta the grid spacing.
        This operator is intended to go be applied to a field with ghosts.
        """
        start = time.time()
        Nc, Nr = 1, 1
        for i in Ns:
            if i > 1:
                Nr *= i-2
                Nc *= i
        Laplacian = csr_matrix((Nr, Nc))
        n = []
        ActualDimension=len([i for i in Ns if i > 1])
        LaplacianDMinusOne=None
        for dir in range(ActualDimension):
            s = dir3D[dir]
            NP2 = list(Ns)
            GL = self._grad(dir, NP2, delta[dir])

            NP1 = list(Ns)
            NP1[dir] += -1
            GR = self._grad(dir, NP1, delta[dir])

            G2 = GR.dot(GL)
            N = list(Ns)
            N[dir] += -2
            ind=self._stripInactiveRows(dir, N)
            Laplacian += G2[ind,:]
            n.append(N[dir])
            if (dir == ActualDimension - 2):
                LaplacianDMinusOne=csr_matrix(Laplacian)


        self._addTime('LaplacianConst', time.time()-start)

        return (Laplacian,LaplacianDMinusOne,n)

    def _restriction(self, RefRatio=None, Ns=None, scaled=True):
        """ define the matrix that provides the restriction operator given a RefRatio """
        start=time.time()
        if RefRatio is None:
            RefRatio=self.RefRatio
        if Ns is None:
            Ns=self._Laplacian[2] # note we use this since we act on interior points always
        if len(Ns)==1:
            return self._restriction((RefRatio[0],1,1), (Ns[0],1,1), scaled)
        if len(Ns)==2:
            return self._restriction((RefRatio[0], RefRatio[1],1), (Ns[0], Ns[1], 1), scaled)
        if len([i for i in RefRatio if i<1])>0:
            raise ValueError("Refratios must be greater than one")


        Nc=1 # accumulates the number of columns
        Nr=1

        for (i,r) in zip(Ns,RefRatio):
            if i>1:
                Nc*=i # number oc columns
                Nr*= i//r # number of rows


        if Nc == Nr: #return an emtpy matrix, since in this case restriction is not used
            A = csr_matrix((1, 1))
            return A,self.GridSize

        if not scaled:
            data=np.ones(Nc)
        else:
            data=np.ones(Nc) # eventually, this must be filled with 1/Jinv

        indicesC=np.arange(Nc)

        indicesR=np.ravel_multi_index((np.unravel_index(indicesC,Ns,order='F')[0]//RefRatio[0],
                                       np.unravel_index(indicesC,Ns,order='F')[1]//RefRatio[1],
                                       np.unravel_index(indicesC,Ns,order='F')[2]//RefRatio[2]),
                                       (Ns[0]//RefRatio[0],
                                       Ns[1]//RefRatio[1],
                                       Ns[2]//RefRatio[2]), order='F')


        A=coo_matrix((data, (indicesR, indicesC)), shape=(Nr,Nc)).tocsr()
        #normalize if scaled
        if scaled:
            ones=np.ones(A.shape[1])
            Dv=A.dot(ones)

            for (i,r) in enumerate(indicesR):
                data[i]/=Dv[r]


            del A
            A=csr_matrix((data,(indicesR, indicesC)), shape=(Nr,Nc))

            self._addTime("Restriction", time.time()-start) # we assume that if not scaled, the operator is
                                                            # called as part of the calculation of prolongation,
                                                            # which times itself.
        #print(" The restriction operator is a ", Nr, " x ", Nc, " matrix.")
        return A,self.GridSize

    def _prolongation(self):
        """ define the prolongation matrix as the transpose of the restriction"""

        start=time.time()
        Ns=[]
        for i,r in zip(self.GridSize, self.RefRatio):
            Ns.append(i*r)
        # note that transpose returns a csc matrix. Untold grief came from this feature not being
        # documented in scipy.
        RV=self._restriction(RefRatio=self.RefRatio, Ns=list(Ns), scaled=False)[0].transpose().tocsr()
        self._addTime("Prolongation", time.time()-start)
        return RV,self.GridSize









if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    import scipy.linalg as la

    Ns = (66, 130, 34)
    Lx, Ly, Lz = (np.pi, 2*(np.pi), np.pi)
    dx, dy, dz = Lx/(Ns[0]-2), Ly/(Ns[1]-2), Lz/(Ns[2]-2)
    DX = (dx, dy, dz)

    x = np.linspace(-dx/2, Lx+dx/2, Ns[0])
    y = np.linspace(-dy/2, Ly+dy/2, Ns[1])
    z = np.linspace(-dz/2, Lz+dz/2, Ns[2])
    yy, xx, zz = np.meshgrid(x, y, z, sparse=False)
    y2d, x2d = np.meshgrid(x, y, sparse=False)

    Phi = np.cos(xx)*np.sin(yy)*np.cos(zz)
    Phi2D = np.cos(x2d) * np.cos(y2d)
    Phi1D = np.cos(x)


    PhiInt = Phi[1:, 1:,:-1]
    V = []
    V.append(Phi[1:, 1:-1, 1:-1])
    V.append(Phi[1:-1, 1:, 1:-1])
    V.append(Phi[1:-1,1:-1,1:])
    Phi1DInt = Phi1D[1:-1]
    Phi2DInt = Phi2D[1:,:-1]


    PhiBox = ((0, 0, 0), Phi.shape, (0, 0, 0))
    Phi2DBox = ((0, 0), Phi2D.shape, (0, 0))
    Phi1DBox = ((0,), Phi1D.shape, (0,))


    NoBC3D = Ellikit(PhiBox, (dy,dx,dz)) # it is dy,dx,dz because of the dumb way meshgrid switches the first two columns...
    NoBC2D = Ellikit(Phi2DBox, (dy, dx))
    NoBC1D = Ellikit(Phi1DBox, (dx, )) # note that the dx's must be in a tuple...




    #plt.plot(x[1:-1], NoBC1D.Dot(Phi1D, 'Laplacian')+Phi1D[1:-1])
    #plt.plot(y[1:-1], NoBC2D.Dot(Phi2D, 'Laplacian')[:, 4]+2*Phi2D[1:-1,5])
    #plt.plot(y[1:-1], NoBC3D.Dot(Phi, 'Laplacian')[:, 4, 8] + 3 * Phi[1:-1, 5, 9])
    print("time in Laplacian without BC", NoBC3D.time['LaplacianConst'])
    PhiBox = ((0, 0, 0), PhiInt.shape, (0, 0, 0))
    Phi2DBox = ((0, 0), Phi2DInt.shape, (0, 0))
    Phi1DBox = ((0,), Phi1DInt.shape, (0,))


    BC3D = Ellikit(PhiBox, (dy,dx,dz), BC=('Neumann', None, 'Dirichlet', None,None, 'Neumann')) # it is dy,dx,dz because of the dumb way meshgrid switches the first two columns...
    print("time in Laplacian with BC", BC3D.time['LaplacianConst'])
    print(" Are the indices sorted ", BC3D.CSRData('Laplacian')['AreIndicesSorted'])
    print(" shape of BC3D grid", BC3D.GridSize)

    BC2D = Ellikit(Phi2DBox, (dy, dx), BC=('Neumann', None, None, 'Neumann'))
    BC1D = Ellikit(Phi1DBox, (dx, ), BC=('Neumann', 'Neumann'), Domain=((0,), Phi1D.shape, (0,))) # note that the dx's must be in a tuple...
    print("offsets ", BC1D.CSRData(Op='Laplacian')['indptr'])
    print('nnz', BC1D.CSRData(Op='Laplacian')['nnz'])

    #plt.plot(x[1:-1], BC1D.Dot(Phi1DInt, 'Laplacian')+Phi1D[1:-1])
    #plt.plot(y[1:-1], BC2D.Dot(Phi2DInt, 'Laplacian')[:, 4] + 2 * Phi2D[1:-1, 5])

    #plt.plot(y[1:-1], BC3D.Dot(PhiInt, 'Laplacian')[:, 4, 8] + 3 * Phi[1:-1, 5, 9])
    print(BC1D._Laplacian)
    print(BC1D.Dot(Phi1DInt, 'Laplacian'))
    print("l2 error 1D ", la.norm(BC1D.Dot(Phi1DInt, 'Laplacian')+Phi1D[1:-1])/Phi1D[1:-1].size)
    print("l2 error 2D ", la.norm(BC2D.Dot(Phi2DInt, 'Laplacian') + 2*Phi2D[1:-1,1:-1])/Phi2D[1:-1,1:-1].size)
    print("l2 error 3D ", la.norm(BC3D.Dot(PhiInt, 'Laplacian') + 3*Phi[1:-1,1:-1,1:-1])/Phi[1:-1,1:-1,1:-1].size)
    print("time in dot Laplacian with BC", BC3D.time['LaplacianDot'])
    print("l2 error 3D of InnerOuter", la.norm(BC3D.Dot(PhiInt, 'Inner')+BC3D.Dot(PhiInt, 'Outer')-PhiInt)/PhiInt.size)

    #plt.plot(x[1:], BC3D.Dot(PhiInt, 'Inner')[2, :, 2])
    #plt.plot(x[1:], BC3D.Dot(PhiInt, 'Outer')[2,:, 2])
    #plt.plot(x[1:], PhiInt[2,:, 2])


    A=BC3D.Dot(PhiInt, 'Inner')+BC3D.Dot(PhiInt, 'Outer')-PhiInt
    plt.show()
    #A=BC3D.Dot(PhiInt,'Outer')
    ind = np.arange(A.size)
    for ix, x in zip(ind, A.ravel('F')):
        if abs(x) > 0:
            k = ix // (A.shape[0] * A.shape[1])
            j = (ix - k * A.shape[0] * A.shape[1]) // A.shape[0]
            i=ix-j*A.shape[0]- k * A.shape[0] * A.shape[1]

            print((i,j,k))



    ErrorNGOG=la.norm(BC3D.Dot(PhiInt[:-1,:-1,1:],'LaplacianNG')+BC3D.Dot(PhiInt,'LaplacianOG')-BC3D.Dot(PhiInt,'Laplacian'))
    print("l2 3D of NG+OG minus total ", ErrorNGOG)


    PhiInt=Phi[1:-1,1:-1,1:-1]
    PhiBox=((0, 0, 0), PhiInt.shape, (0, 0, 0))

    div=BC3D.Dot(V,'Divergence')

    BC3D=Ellikit(PhiBox, (dy,dx,dz), RefRatio=(2,2,2))
    #plt.plot(x[1:-1:2], BC3D.Dot(PhiInt, Op="Restriction")[1,:,1])
    #plt.plot(x[1:-1], PhiInt[2,:,2])
    #plt.show()

    #print(NoBC3D._MainDiagonals[0])
    #print("---")
    #print(NoBC3D._MainDiagonals[1])

    print(" ***** Beginning test of reconstructed operator *********************")

    # these are the boxes overlapping on the interior points.
    BC0 = ('Neumann', None, 'Neumann', 'Neumann')
    BC1 = (None, 'Neumann', 'Neumann', None)
    BC2 = ('Neumann', 'Neumann', None, 'Neumann')
    BC = (BC0, BC1, BC2)
    Phi2D=[]
    Box = []
    D=[]
    RefRatio = ((2,2),(1,1))
    A2D = []
    LDV = []
    S=[]
    nx, ny = 32, 16 # size of individual boxes

    for lev,RR in enumerate(RefRatio):
        B0 = ((0, 0), (nx+1, ny), (0, 0))

        B1 = ((nx-1, 0), (nx+1, ny+1), (0, 0))

        B2 = ((nx, ny-1), (nx, ny+1), (0, 0))

        D.append(((0, 0), (2 * nx, 2 * ny), (0, 0)))  # this is the enveloping domain.
        Delta=(2*np.pi/D[lev][1][0], 2*np.pi/D[lev][1][1])
        dx,dy=Delta
        x = np.linspace(dx/2, 2*np.pi-dx/2, D[lev][1][0])
        y = np.linspace(dy/2, 2*np.pi-dy/2, D[lev][1][1])

        x2d, y2d = np.meshgrid(y, x, sparse=False)


        Phi2D.append(np.cos(x2d) * np.cos(y2d)) # an eigenvalue of the operator

        print('Phi2D', Phi2D[lev].shape)
        Box.append((B0, B1, B2))


        A2D.append([])
        for b, bc in zip(Box[lev], BC):

            A2D[lev].append(Ellikit(b, Delta, BC=bc, Domain=D[lev], RefRatio=RR))



    # accumulate the full laplacian. Note that at this point the container is defined over
    # the enveloping domain.
        L = csr_matrix((4*nx*ny,4*ny*nx))
        for A in A2D[lev]:
            L += A._G2D[0].dot(A._Laplacian[0].dot(A._D2G))




    #plt.imshow(L.dot(Phi2D.ravel('F')).reshape(D[1], order='F'))

    #plt.show()


    # pack a LD equivalent

        LDV.append((Phi2D[lev][:nx,:ny], Phi2D[lev][nx:,:ny], Phi2D[lev][nx:,ny:]))
        S.append(csr_matrix((L.shape[0],1))) # this is the container over the enveloping domain.
        for V, A in zip(LDV[lev], A2D[lev]):
        # we pack data over a box into a (n X 1) sparse matrix so that grid to domain returns
        # a sparse matrix (instead of the full vector)

            v = csr_matrix((V.ravel('F'), np.zeros(nx * ny, dtype='int'), np.arange(nx*ny + 1)))

            S[lev] += A._G2D[0].dot(v)

        #S.indptr need to be saved as it will be used later to unpack the full solution
        # over individual grids.
        plt.imshow(S[lev].A.reshape(D[lev][1], order='F'))
        plt.colorbar()
        plt.show()

    # From the enveloping domain, we extract the field over the actual domain.
        mask=np.concatenate(([True], S[lev].indptr[1:]!=S[lev].indptr[:-1]))
        SNonZero = csr_matrix((S[lev].data, S[lev].indices, S[lev].indptr[mask]))

    # restrict the Laplacian to the actual domain,
    #eliminating the nonzero rows. Note we use the same mask derived from S
        #mask = np.concatenate(([True], L.indptr[1:] != L.indptr[:-1]))
        L = csr_matrix((L.data, L.indices, L.indptr[mask]))
        print("shape of L after first pass ", L.shape)
    # now we elimiate the nonzero columns.
    # We do this transposing the matrix and working on the rows of the transposed matrix.
    # Note that scipy.sparse.csr_matrix.transpose()
    # returns a csc_matrix. Hence the need for the conversion.
        L = L.transpose().tocsr()
        mask = np.concatenate(([True], L.indptr[1:] != L.indptr[:-1]))
        L = csr_matrix((L.data, L.indices, L.indptr[mask]))
    # this step is likely not necessary, as L at this point should be symmetric...
        L=L.transpose().tocsr()




    # solve equation L to reduced vector
    # notice that here SNonZero is a n X 1 sparse matrix.
    #LS = L.dot(SNonZero)

        LS,info=LA.cg(L,-2*SNonZero.A,tol=1.e-5,  maxiter=10000)
    # expand  LS onto the full grid. We use the indptr saved from the previous ste.
        LFull = csr_matrix((LS, np.zeros(LS.shape[0], dtype='int'), S[lev].indptr))
        # note that we use .A to extract all components from the nX1 matrix, even those that may be zero.
        plt.imshow(LFull.A.reshape(D[lev][1], order='F'))
        plt.colorbar()
        plt.show()
    # repack the data into a LD

        LofLDV=3*[None]
        for n, (A, b) in enumerate(zip(A2D[lev], Box[lev])):

            V = A._D2G.dot(LFull)

            LofLDV[n] = V.A.reshape(b[1], order='F')

        for n, (Lv, v) in enumerate(zip(LofLDV, LDV[lev])):

            plt.imshow(Lv[:nx,:ny])
            plt.colorbar()
            plt.show()

        nx, ny = nx // RR[0], ny // RR[1]

    Restriction = len(RefRatio)*[None]
    Prolongation=len(RefRatio)*[None]
    for lev in range(len(RefRatio)):
        if (lev < len(RefRatio)-1):
            # work on the restriction operator
            R = csr_matrix((D[lev + 1][1][0] * D[lev + 1][1][1], D[lev][1][0] * D[lev][1][1]))
            for AFine, ACoarse in zip(A2D[lev], A2D[lev + 1]):

                R += ACoarse._G2D[0].dot(AFine._Restriction[0].dot(AFine._G2D[0].transpose()))
            mask=np.concatenate(([True], S[lev+1].indptr[1:]!=S[lev+1].indptr[:-1]))
            Restriction[lev]=R #csr_matrix((R.data, R.indices, S[lev + 1].indptr[mask]))

    for lev in range(len(RefRatio)):

        if (lev > 0):

            Prolongation[lev] = RefRatio[lev-1][0]*RefRatio[lev-1][1]*Restriction[lev-1].transpose().tocsr()



    Restricted = Restriction[0].dot(S[0])

    plt.imshow(Restricted.A.reshape(D[1][1], order='F'))
    plt.colorbar()
    plt.show()

    Prolongated = Prolongation[1].dot(S[1])
    plt.imshow(Prolongated.A.reshape(D[0][1], order='F'))
    plt.colorbar()
    plt.show()

