#!/usr/bin/env python3
""" generates a LevelDataFAB analytically."""
from LevelDataFAB import LevelDataFAB as ldFAB
from LevelDataFB import LevelDataFluxBox as LDFB
from FAB import FAB
from Box import Box as Box
import Slicer
import Diagnostic
from random import seed
seed(1)
from random import randint
import matplotlib.pyplot as plt
from mpi4py import MPI
import numpy as np
def dimension():
    return 3
def domainLength():

    return (.2,.3,.4)
def domainSize():
    return (32,64,128)
def minBlock():
    #return (domainSize())
    return (16,16,64)

def myFunction(LDCoord,LDField):
    for x,f in zip(LDCoord,LDField):
        for d in range(3):
            
            X=x.Fluxes[d][...,0]
            Y=x.Fluxes[d][...,1]
            Z=x.Fluxes[d][...,2]
            if d==0:f.Fluxes[d][...,0]=X*X+Y*Y
            if d==1:f.Fluxes[d][...,0]=Y*Y+Z*Z
            if d==2:f.Fluxes[d][...,0]=Z*Z+X*X
            #f.Fluxes[d][...,0]=np.tanh(2*(Z-6.95/2))/2
        
        
    return 
def myFunctionFAB(LDCoord, LDField):
    for x,f in zip(LDCoord,LDField):
        X=x[...,0]
        Y=x[...,1]
        Z=x[...,2]
        f[...,0]=X*Y*Z
        f[...,1]=0.
        f[...,2]=0.
    return


if __name__ == "__main__":
   
    # seed random number generator
    
    
    n=domainSize() # total domain
    m=minBlock()  # min block
    q=tuple([S//B for S,B in zip(n,m)])
    NBoxes=q[0]
    factors=dimension()*[1]
    for d,(c,cprevious) in enumerate(zip(q[1:],q[0:-1])):
        NBoxes*=c
        factors[d+1]=factors[d]*cprevious

   
    L=domainLength() # domain physical size
    # create collection of boxes
    comm = MPI.COMM_WORLD
    MyRank = comm.Get_rank()
    MySize = comm.Get_size()

    Boxes = []
    Pids = []
    FBs = []
    for c in range(NBoxes):
        LoEnd=(dimension()+1)*[0]
        r=c
        
        for d in range(dimension()-1,-1,-1):
            k=r//factors[d]
            r-=k*factors[d]
            LoEnd[d]=k*m[d]
        if r!=0: raise(ValueError("Remainder not zero"))
        
        b = Box.FromLoEndAndSize(tuple(LoEnd[:-1]), m)
        
        Boxes.append(b)
        Pids.append((randint(0,MySize-1),))

    # Generate the boxes that cover the domain and assign them to a random process
   # for k in range(n[2] // m[2]):
   #     for j in range(n[1] // m[1]):
   #         for i in range(n[0] // m[0]):
   #             b = Box.FromLoEndAndSize((i * m[0], j * m[1], k * m[2]), m)
   #             Boxes.append(b)

   #             Pids.append((randint(0, MySize - 1), ))


    # create fbs
    
    
    for b, id in zip(Boxes, Pids):

        if id[0] == MyRank:
            Args = []

            G = []
            S = []
            for l, s in zip(b.LoEnd(), b.Size()):
                G.append(l - 1)
                S.append(s + 2)

            Args.append(tuple(G)) # low end
            Args.append(tuple(S)) # size
            Args.append(dimension()*(0,)) # centering
            Args.append((dimension(),)) # 3 components
            Args.extend(dimension()*[None]) # allocate memory
            Args.append('FluxBox')
            FBs.append(Args)

    Args=5*[None]
    Args[4] = 'LevelDataFluxBox'

    Args[3] = dimension()*(0,) +tuple([s-1 for s in n]) + ('Box',)
    Args[2] = Pids
    Args[1] = [(b.LoEnd() + b.HiEnd()+('Box',)) for b in Boxes]
    Args[0] = FBs
    LDCoord = LDFB(Args)

    # create fbs
    Args=[]
    FBs=[]
    for b, id in zip(Boxes, Pids):

        if id[0] == MyRank:
            Args = []

            G = []
            S = []
            for l, s in zip(b.LoEnd(), b.Size()):
                G.append(l - 1)
                S.append(s + 2)

            Args.append(tuple(G)) # low end
            Args.append(tuple(S)) # size
            Args.append(dimension()*(0,)) # centering
            Args.append((1,)) # 3 components
            Args.extend(dimension()*[None]) # allocate memory
            Args.append('FluxBox')
            FBs.append(Args)

    Args=5*[None]
    Args[4] = 'LevelDataFluxBox'

    Args[3] = dimension()*(0,) +tuple([s-1 for s in n]) + ('Box',)
    Args[2] = Pids
    Args[1] = [(b.LoEnd() + b.HiEnd()+('Box',)) for b in Boxes]
    Args[0] = FBs

    LDField = LDFB(Args)
    if MyRank==0: print("done creating the LDFBs")
    # create a LDFAB on the same box layout
    # create a fabs with ghost layers
    FABs=[]
    for b, id in zip(Boxes, Pids):

        if id[0] == MyRank:
            Args = []

            G = []
            S = []
            for l, s in zip(b.LoEnd(), b.Size()):
                G.append(l - 1)
                S.append(s + 2)

            Args.append(tuple(G)) # low end
            Args.append(tuple(S)) # size
            Args.append(dimension()*(0,)) # centering
            Args.append((dimension(),)) # 3 components
            Args.append(None) # allocate memory
            Args.append('FAB')
            FABs.append(Args)

    Args=5*[None]
    Args[4] = 'LevelDataFAB'

    Args[3] = dimension()*(0,) +tuple([s-1 for s in n]) + ('Box',)
    Args[2] = Pids
    Args[1] = [(b.LoEnd() + b.HiEnd()+('Box',)) for b in Boxes]
    Args[0] = FABs

    Domain=Box(Args[3])
    LDFABField = ldFAB(Args)
    
    LDFABCoord = ldFAB.MakeFromLDFAB(LDFABField)
    

    
    perm=((0,1,2),(1,2,0),(2,0,1)) if dimension()==3 else ((0,1),(1,0))
    

    #fill the coordinates [0,Lx]x[0,Ly]x[0,Lz]
    integral=0
    delta=[l/s for l,s in zip(L,n)]
    volume=1
    for dx in delta:
        volume*=dx
    for x in LDCoord:
        
        
        for d,f in enumerate(x.Fluxes):
            b=f.box
            #print(b.Centering())
            for p in perm:
                start=(b.LoEnd()[p[0]])
                end=start+b.Size()[p[0]]
            
            
                replicate=1
                for s in p[1:]:
                    replicate*=b.Size()[s]
                newshape=tuple([b.Size()[s] for s in p])
                offset=0.5 if f.centering[p[0]]==0 else 0
                f[...,p[0]]=np.moveaxis(np.tile((np.arange(start,end)+offset)*delta[p[0]],replicate).reshape(newshape,order='F'),perm[0],p)   

            #if (d==0): print(f[:,1,1,0])
            #if (d==1): print(f[1,:,1,1])
            #if (d==2): print(f[:,1,1,0])
   
    #fill the LDFAB
    #     
    myFunction(LDCoord,LDField) # fills the field
    if MyRank==0: print("Initialized LDField")
    integral=0
    delta=[l/s for l,s in zip(L,n)]
    volume=1
    for dx in delta:
        volume*=dx
    for x in LDFABCoord:
        
        b=x.box
        for p in perm:
            start=(b.LoEnd()[p[0]])
            end=start+m[p[0]]+2
            
            
            replicate=1
            for s in p[1:]:
                replicate*=(m[s]+2)
            newshape=tuple([m[s]+2 for s in p])
            
            x[...,p[0]]=np.moveaxis(np.tile((np.arange(start,end)+0.5)*delta[p[0]],replicate).reshape(newshape,order='F'),perm[0],p)   


    myFunctionFAB(LDFABCoord,LDFABField) # fills the field
    if MyRank==0: print("Initialized LDFABField")
    #acc=0
    #for i in range(3):
    #    for j in range(i,3):
    #        intd11=0
    #        for f in LDField:
    #            d11=(Diagnostic.differentiateFAB(f,delta,i,j)+Diagnostic.differentiateFAB(f,delta,j,i))/2
    #            intd11+=2*np.sum(d11*d11)*volume if (i==j) else 4*np.sum(d11*d11)*volume

    #        s=MPI.COMM_WORLD.allreduce(intd11,op=MPI.SUM)
    #        print(">>>",i,j,s)
    #        acc+=s
    #print(MPI.COMM_WORLD.Get_rank(),acc)
    #print(delta)
    nu=1.5
    KEAnalytic=(L[0]*L[1]*L[2])*((L[0]**4+L[1]**4+L[2]**4)/5+(L[0]**2*L[1]**2+L[1]**2*L[2]**2+L[2]**2*L[0]**2)/9)
    
    # EpsilonAnalytic=nu*(4*L[0]*L[1]*L[2]*(L[0]**2+L[1]**2+L[2]**2))
    EpsilonAnalytic=.75*(L[0]**2+L[1]**2+L[2]**2)*(L[0]*L[1]*L[2])**2
    DeltaKE=abs(Diagnostic.KE(LDField,delta)-KEAnalytic)/KEAnalytic
    if MyRank==0: print("Done with KE")
    Epsilon,LocalEpsilon,_=Diagnostic.epsilonSGS(LDField,LDFABField,delta)
    if MyRank==0: print("Done with epsilon")
    if(MPI.COMM_WORLD.Get_rank()>=0):
        print(Epsilon,LocalEpsilon)
        print(abs(Epsilon-EpsilonAnalytic)/Epsilon)
        
        

    
    




        