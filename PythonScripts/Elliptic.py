

import PyGlue as pg
#from numba import jit
from mpi4py import MPI
from Box import Box
from FluxBox import FluxBox
from LevelDataFB import LevelDataFluxBox
from LevelDataFAB import LevelDataFAB
from SOMAR import SOMAR
from numba import jit
import numpy as np
@jit(nopython=True)
def execute(phi, res, invdiags, omega, whichPass, iRange,
                                                 phiOff,
                                                 resOff,
                                                 invdiagsOff):

    for j in range(iRange[2], iRange[3]+1):

        indtot=iRange[0]+j
        start=iRange[0]+abs((indtot+whichPass)%2)

        for i in range(start, iRange[1]+1,2):

            phi[i-phiOff[0],j-phiOff[1]]+=omega*invdiags[i-invdiagsOff[0],j-invdiagsOff[1]]*res[i-resOff[0],j-resOff[1]]

@SOMAR
def PoissonOp_JacobiRB(phi, res, invDiags, omega, whichPass, destBox):
    if phi.box.ContainsBox(destBox):
        execute(phi.data[...,0], res.data[...,0], invDiags.data[...,0], omega, whichPass,
                np.array([destBox.LoEnd()[0],destBox.HiEnd()[0],
                destBox.LoEnd()[1], destBox.HiEnd()[1]]),
                np.array(phi.offsets),
                np.array(res.offsets),
                np.array(invDiags.offsets))
    else:
        raise ValueError





