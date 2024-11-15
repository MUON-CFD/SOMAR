import Slicer
import os
from SOMAR import SOMAR
import Comm
import string
import random
from mpi4py import MPI

def _id_generator():
    chars=string.ascii_letters + string.digits
    return str().join(random.SystemRandom().choice(chars) for _ in range(10))

@SOMAR
def PlotPlane(LD, fileName, title, dir, pos, comp):
    
    Slicer.PlotSlice(LD, dir= dir, pos = pos, comp = comp, filename = fileName, title = title)


@SOMAR
def PlotLine(LD, fileName, title, dir, pos, comp):
    Slicer.Plotline(LD, dir = dir, pos = pos, comp = comp, filename = fileName, title = title)

@SOMAR 
def MailPlotPlane(LD, address, title, dir, pos, comp):
    
    fileName = _id_generator()+'.pdf'
    Slicer.PlotSlice(LD, dir= dir, pos = pos, comp = comp, filename = fileName, title = title)
    if LD.myMpiId==0:
        Comm.sendMailTo(address, filename = fileName)
        os.remove(fileName)


@SOMAR
def MailPlotLine(LD, address, title, dir, pos, comp):

    fileName = _id_generator()+'.pdf'
    Slicer.Plotline(LD, dir = dir, pos = tuple(p for p in pos), comp = comp, filename = fileName, title = title)
    if LD.myMpiId==0:
        Comm.sendMailTo(address, filename = fileName)
        os.remove(fileName)

    