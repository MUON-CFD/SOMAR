""" a barebone PyPhysics class, to be used in case we make CH_USE_PYTHON the default,
 but we do not need Python services. """


import os
from os import path
import sys
CWD = os.getcwd()
PyScriptsDir = CWD + '/../PythonScripts'

sys.path.append(PyScriptsDir)
import numpy as np
import re
try:
    import matplotlib.pyplot as plt
except ImportError:
    print(" Matplotlib not available. Consider installing it.")
import math
import PyGlue as pg
#from numba import jit
from mpi4py import MPI
from FluxBox import FluxBox
from LevelDataFB import LevelDataFluxBox
from LevelDataFAB import LevelDataFAB
from SOMAR import SOMAR



class Physics:
    """ stores value common to all processes. """
    parmDict = dict()
    compDict = dict()
    is_Initialized = False


    comm = MPI.COMM_WORLD
    MpiNumProcs = comm.Get_size()

    def __init__(self): pass



    @classmethod
    def GetNumProcs(cls):
        
        return cls.MpiNumProcs

    @classmethod
    def GetComm(cls):
        
        return cls.comm

    @classmethod
    def GetParameters(cls):
        
        return cls.parmDict

    @staticmethod
    def strType(var):
        try:
            if int(var) == float(var):
                return int(var)
        except:
            try:
                float(var)
                return float(var)
            except:
                return var

    @classmethod
    def process_input(cls,ss):
        
        for i in range(len(ss)):
            for ch in ['(', ')', '=', ',']:
                if ch in ss[i]:
                    ss[i] = ss[i].replace(ch, '')

        tup = ()
        t_l = []
        i = 0
        while i < len(ss):
            num = int(re.findall(r'\d+', ss[i])[0])
            if num == 1:
                val = cls.strType(ss[i+3])
                tup = (ss[i+1], val)
            else:
                temp = []
                for k in range(num):
                    val = cls.strType(ss[i+3+k])
                    temp.append(val)
                tup = (ss[i+1], temp)
            t_l.append(tup)
            i = i + num + 3
        for x, y in t_l:
            cls.parmDict[x] = y

    @classmethod
    def ReadParms(cls,inputfile='inputs.template.machine'):

        try:
            f=open(inputfile,'r')
            parms=f.read()
            f.close()
            ss=parms.split('\n') # split lines
            x=[s.split() for s in ss if (len(s)>0) and s[0]!='#'] # eliminate comment lines
            for line in x:
                y=[]
                j=[i for i,entry in enumerate(line) if entry=='#'] # only consider entries up to a #
                if len(j)>0:
                    ll=line[2:j[0]]
                else:
                    ll=line[2:]
                for entry in ll:
                    try:
                        y.append(int(entry))
                    except:
                        try:
                            y.append(float(entry))
                        except:
                            y.append(entry)
                    finally:
                        cls.parmDict[line[0]] = y[0] if len(y)==1 else tuple(y)


        except:
            raise ValueError("Cannot open input file named %s " % inputname)

    @classmethod
    def Initialize(cls,inputName, MpiId):  # scalar names  are defined above
        
        cls.myMpiId = MpiId

        if cls.myMpiId != cls.comm.Get_rank():
            raise ValueError("The MPI Rank passed by C++ does not correspond to what I see here.")


        if cls.is_Initialized: return
        cls.is_Initialized = True
        cls.ReadParms(inputName)



@SOMAR
def Initialize(inputName, MyMpiId):
    Physics().Initialize(inputName, MyMpiId)



@SOMAR
def GetParm(name):

    try:
        return Physics().parmDict[name]
    except:
        raise ValueError("Parameter %s is not defined " % name)
@SOMAR
def Conclude(): pass
