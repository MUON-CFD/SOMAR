import numpy as np
from FAB import FAB as FAB
from LevelDataFAB import LevelDataFAB as LDFAB
from LevelDataFB import LevelDataFluxBox as LDFB
from mpi4py import MPI
from Box import Box
from SOMAR import SOMAR
import UserDefinedBoundaryFunctions as FillFunc

VDiri=FillFunc.DiriVelocity()
VNeumann=FillFunc.NeumannVelocity()

@SOMAR
def SetVelAtBoundary(BufferToBeFilled, VelDirection, BoxOverWhichWeFill, BdryDir, BdrySide, CoordsOverBox, time):
    VDiri.Fill(BufferToBeFilled, VelDirection, BoxOverWhichWeFill, BdryDir, BdrySide, CoordsOverBox, time)


