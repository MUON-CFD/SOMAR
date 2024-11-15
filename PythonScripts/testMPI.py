from mpi4py import MPI

print("My rank is ", MPI.COMM_WORLD.Get_rank())

print("My size is ", MPI.COMM_WORLD.Get_size())
import h5py as h5
f=h5.File('test.h5','a',driver='mpio',comm=MPI.COMM_WORLD)

f.close()
