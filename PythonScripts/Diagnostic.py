"""Sundry collection of diagnostic utilities"""
from FluxBox import FluxBox as FluxBox
from FAB import FAB as FAB
from SOMAR import SOMAR
import numpy as np
from mpi4py import MPI
#from multipledispatch import dispatch

@SOMAR
def stats(LDFAB):
    M=0.
    ave=0.
    Ncells=0.
    Slice=tuple(slice(1,-1) for d in range(LDFAB.SpaceDim))
    for box in LDFAB.Boxes:
        Ncells+=box.Ncells()

    for fab in LDFAB:
        valid=fab[Slice]
        
        M=max(abs(valid[...]).max(),M)
        ave+=np.sum(valid[...])/Ncells
        
    M=MPI.COMM_WORLD.allreduce(M,op=MPI.MAX)
    ave=MPI.COMM_WORLD.allreduce(ave,op=MPI.SUM)
    std=0.
    for fab in LDFAB:
        valid=fab[Slice]
        
        std+=np.sum(np.square((valid[...]-ave)))/Ncells
        
    std=MPI.COMM_WORLD.allreduce(std,op=MPI.SUM)
    return (M,ave,np.sqrt(std))

@SOMAR
def Richardson(V, B, RI, COORD):
    """Given a Velocity Field V, a Buoyancy field B, returns the Richardson number Ri. 
    All inputs are assumed to be LevelDatas. To make is so that it works in multiple
    dimension, dir must be specified  """
    epsilon = 1.0e-12
    for (ri, v, b, coord) in zip(RI, V, B, COORD):
        SpaceDim = len(b[..., 0].shape)
        u = (v.Fluxes[0][1:-2, ...]+v.Fluxes[0][:-3, ...])/2
        z = coord[..., SpaceDim-1]

        if SpaceDim == 3:
            dz = z[1, 1, 2]-z[1, 1, 0]
            S2 = ((u[:, 1:-1, 2:] - u[:, 1:-1, 0:-2]) / (dz)) ** 2
            N2 = 1-((b[1:-1, 1:-1, 2:, 0] - b[1:-1, 1:-1, 0:-2, 0]) / (dz))

        else:
            dz = z[1, 2]-z[1, 0]
            S2 = ((u[:, 2:] - u[:, 0:-2]) / (dz)) ** 2
            N2 = 1-((b[1:-1, 2:, 0] - b[1:-1, 0:-2, 0]) / (dz))

        ri[..., 0] = N2 / (S2[..., 0]+epsilon)

@SOMAR
def PE(B,DX,Z,bComp=0):
    localPE=0.
    SpaceDim=B.SpaceDim
    dx=3*[1.]
    dx[0:SpaceDim]=DX[0:SpaceDim]
    dV=1.
    Slice=tuple(slice(1,-1) for d in range(SpaceDim))
    for d in dx:
        dV*=d
    for b,z in zip(B,Z):
        PE=b[Slice+(bComp,)]*z[Slice+(SpaceDim-1,)] 
        localPE+=np.sum(PE)*dV 
    
    return MPI.COMM_WORLD.allreduce(localPE,op=MPI.SUM)

@SOMAR
def BF(V,Q,DX, bcomp=0):
    D = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    DC= ((0, 1, 1), (1, 0, 1), (1, 1, 0))
    LocalBF=0.
    SpaceDim=V.SpaceDim
    n=3*[1]
    dx=3*[1.]
    dx[0:SpaceDim]=DX[0:SpaceDim]
    dV=1.
    for d in dx:
        dV*=d
    Slice=tuple(slice(1,-1) for d in range(SpaceDim))+(0,)
    for v,q in zip(V,Q):
        w=v.Fluxes[SpaceDim-1][...]
        n[0:v.dimension]=w.shape[0:v.dimension]
        # bring the FC to CC 
        
        w=(w[      D[SpaceDim-1][0]:,
                   D[SpaceDim-1][1]:,
                   D[SpaceDim-1][2]:]+
           w[:n[0]-D[SpaceDim-1][0],
             :n[1]-D[SpaceDim-1][1],
             :n[2]-D[SpaceDim-1][2]])/2
        
        WB=w[Slice]*q[Slice]
        
        LocalBF+=np.sum(WB)*dV
    
    return MPI.COMM_WORLD.allreduce(LocalBF, op=MPI.SUM)

@SOMAR
def BFFAB(Q,DX,wcomp=1,bcomp=5):
    LocalBF=0.
    dx=3*[1]
    dx[0:Q.SpaceDim]=DX[0:Q.SpaceDim]
    dV=1.
    Slice=tuple(slice(1,-1) for d in range(Q.SpaceDim))
    for d in dx:
        dV*=d
    for q in Q:
        w=q[...,wcomp]
        b=q[...,bcomp]
        
        WB=w[Slice]*b[Slice]
        LocalBF+=np.sum(WB)*dV
    return MPI.COMM_WORLD.allreduce(LocalBF, op=MPI.SUM)

def Strain(v,dx):
    # only the upper diagonal part is stored. 
    SpaceDim=v.dimension
    S=SpaceDim*[[]]
    for i in range(SpaceDim):
        R=SpaceDim*[[]]
        R[i]=differentiateFB(v,dx,i,i) if isinstance(v, FluxBox) else differentiateFAB(v,dx,i,i)
              
        for j in range(i+1, SpaceDim):
            R[j] = (differentiateFB(v, dx, i, j)+ \
                               differentiateFB(v, dx, j, i))/2 if isinstance(v, FluxBox) else \
                          (differentiateFAB(v, dx, i, j)+ \
                               differentiateFAB(v, dx, j, i))/2        
                        
                
        S[i]=R
    return S

@SOMAR
def epsilonField(V,nut,out,DX,nu):
    SpaceDim=V.SpaceDim
    dx=[1,1,1]
    dx[0:SpaceDim]=DX[0:SpaceDim]
    Slice=tuple(slice(1,-1) for dir in range(SpaceDim))

    for _v,_nut,_out in zip(V,nut,out):
        S=Strain(_v,dx)
        
        _out[Slice+(0,)]=2 * traceStrainSquaredField(S) * (nu + _nut[Slice+(0,)])

    return

@SOMAR
def ChiField(S,nut, out, DX, kappa):
    SpaceDim = S.SpaceDim
    dx=[1,1,1]
    dx[0:SpaceDim]=DX[0:SpaceDim]
    Slice=tuple(slice(1,-1) for dir in range(SpaceDim))

    for _s,_nut,_out in zip(S,nut,out):
        dummy=9.81 * 8e-4 * differentiateFAB(_s,DX,0,0)
        _out[Slice+(0,)] = (kappa + _nut[Slice+(0,)]) * dummy * dummy
        for d in range(1,SpaceDim):
            dummy=9.81 * 8e-4 * differentiateFAB(_s,DX,0,d)
            _out[Slice+(0,)] += (kappa + _nut[Slice+(0,)]) * dummy * dummy
    
    return


        
@SOMAR 
def epsilonSGS(V,NUT,DX,nuTComp=0):
    """Given a velocity field V (A leveldata) returns the dissipation field 
    \epsilon=\int 2*nu*Tr (S*S), where S is the rate of strain tensor. Assumes cartesian grid
    with spacing dx"""
    LocalTrace = 0.
    LocalDissipation = 0.
    SpaceDim=V.SpaceDim
    dx=[1,1,1]
    dx[0:SpaceDim]=DX[0:SpaceDim]
    dV=1.
    for d in range(SpaceDim):
        dV*=dx[d]
    M=0.
    Slice=tuple(slice(1,-1) for dir in range(SpaceDim))+(nuTComp,)
    for  v,nut in zip(V,NUT): # dit loop
        
        S=Strain(v,dx)
        try:
            eps=traceStrainSquaredField(S)
        except:
            eps=traceStrainSquared(S)
        LocalTrace+= np.sum(2*eps)*dV
        
        nu=nut[Slice]
        
            
        LocalDissipation += np.sum(2*nu*eps)*dV
        

    # finally, accumulate
    totalLocalDissipation = MPI.COMM_WORLD.allreduce(LocalDissipation, op=MPI.SUM)
        
    return (totalLocalDissipation, LocalDissipation, LocalTrace)
        
@SOMAR
def epsilon(V, nu, DX):
    """Given a velocity field V (A leveldata) returns the dissipation field 
    \epsilon=\int 2*nu*Tr (S*S), where S is the rate of strain tensor. Assumes cartesian grid
    with spacing dx"""
    
    LocalDissipation = 0.
    SpaceDim=V.SpaceDim
    dx=[1,1,1]
    dx[0:SpaceDim]=DX[0:SpaceDim]
    dV=1.
    for d in range(SpaceDim):
        dV*=dx[d]

    for  v in  V: # dit loop
        S=Strain(v,dx)
        try:
            eps=traceStrainSquaredField(S)
        except:
            eps=traceStrainSquared(S)
                
        LocalDissipation += np.sum(2*nu*eps)*dV

        

    # finally, accumulate
    totalLocalDissipation = MPI.COMM_WORLD.allreduce(LocalDissipation, op=MPI.SUM)
        
    return totalLocalDissipation
        
def traceStrainSquaredField(S):
    # the components of S are defined on different cell centering
    # we bring everything to the same cell centering after we squared the
    # components 
    # die Kunst der Teil
    D = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    DC= ((0, 1, 1), (1, 0, 1), (1, 1, 0))
    SpaceDim=len(S[0])
    n=[]
    for d in range(SpaceDim):
        n.append(S[0][0].shape[d])
        n[d]-= 2
    
    eps=np.zeros(n)
    
    for i in range(len(S)):
        Slice=tuple(slice(1-D[i][d],-1-D[i][d]) for d in range(SpaceDim))
        
        eps += np.square(S[i][i][Slice])
        
        
        for j in range(i+1,len(S)):
            Slice=[slice(0,None) for d in range(SpaceDim)]
            nx=S[i][j].shape
            
            for m in range(2):
                
                Slice[j]=slice(m,nx[j]-(1-m))
            
                for n in range(2):
                    
                    Slice[i]=slice(n,nx[i]-(1-n))
                    
                    eps += np.square(S[i][j][tuple(Slice)])/2



            
    return np.squeeze(eps)



def traceStrainSquared(S):
    eps=np.zeros_like(S[0][0])
    for i in range(len(S)):
        eps += np.square(S[i][i])
        for j in range(i+1,len(S)):
            eps+=2*np.square(S[i][j])
    return np.squeeze(eps)


def differentiateFB(v , dx, i, j):
    D = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    DC= ((0, 1, 1), (1, 0, 1), (1, 1, 0))
    n=[1,1,1]
    u=v.Fluxes[i][...,0] if v.dimension==3 else v.Fluxes[i][...]
    n[0:v.dimension]=u.shape[0:v.dimension]
    sliceLeft=tuple(slice(D[j][d],None) for d in range(v.dimension))
    sliceRight=tuple(slice(None, n[d]-D[j][d]) for d in range(v.dimension))
    derivative=(u[sliceLeft]-u[sliceRight])/dx[j]
                    
    if(i!=j):
        Slice=tuple(slice(0,None) if d == j else slice(1,-1) for d in range(v.dimension))
        derivative=derivative[Slice]
    
    return derivative if v.dimension==3 else derivative[...,0]


    
def differentiateFAB(v,dx,i,j): 
    D = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    DC= ((0, 1, 1), (1, 0, 1), (1, 1, 0))
    n=[1,1,1]
    n[0:v.dimension]=v.size
    u=np.empty(tuple(n),dtype='double', order='F') 
    stencil=[-1/(4*dx[j]),6/(4*dx[j]),-5/(4*dx[j])]
    #stencil = [1/(2*dx[j]), 0, -1/(2*dx[j])]
    sliceLeft=tuple(slice(2*D[j][d], None) for d in range(v.dimension))
    sliceCenter=tuple(slice(D[j][d],n[d]-D[j][d]) for d in range(v.dimension))
    sliceRight=tuple(slice(None, n[d]-2*D[j][d]) for d in range(v.dimension))
    try:
        u[...]=v[...,i]
    except:
        u[...,0]=v[...,i]
    derivative=(stencil[0]*u[sliceLeft]+stencil[1]*u[sliceCenter]+stencil[2]*u[sliceRight])
    sliceOut = tuple(slice(DC[j][d], n[d]-DC[j][d]) for d in range(v.dimension))
    if v.dimension==2:
        sliceOut+=(0,)
    return derivative[sliceOut]
    

@SOMAR
def KE(V,dx):
    D = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    DC= ((0, 1, 1), (1, 0, 1), (1, 1, 0))
    LocalKE=0.
    dV=1.
    n=[1,1,1]
    Slice=tuple(slice(1,-1) for d in range(V.SpaceDim))+(0,)
    for d in dx:
         if (isinstance(d,(int,float))):
            dV*=d
    
    for q in V:
        for i in range(V.SpaceDim):
            
            u=q.Fluxes[i][...]
            n[0:q.dimension]=u.shape[0:q.dimension]
             # bring the FC to CC 
            uCC=(u[      D[i][0]:,
                   D[i][1]:,
                   D[i][2]:]+
                u[:n[0]-D[i][0],
                    :n[1]-D[i][1],
                :n[2]-D[i][2]])/2
            
            LocalKE+=np.sum(uCC[Slice]*uCC[Slice])/2*dV 
            
    return MPI.COMM_WORLD.allreduce(LocalKE,op=MPI.SUM)

@SOMAR
def KEFAB(V,dx):
    LocalKE=0.
    dV=1.
    Slice=tuple(slice(1,-1) for d in range(V.SpaceDim))+(slice(0,V.SpaceDim),)
    for d in dx:
        if (isinstance(d,(int,float))): dV*=d
    for q in V:
        LocalKE+=np.sum(q[Slice]*q[Slice])/2*dV
        
    return MPI.COMM_WORLD.allreduce(LocalKE,op=MPI.SUM)

