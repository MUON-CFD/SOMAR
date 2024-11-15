''' Provides a decorator to convert the arguments of a function from lists representing objects as passed by the SOMAR calling
function to the corresponding objects in Python.Passes the converted object to the function '''

from LevelDataFAB import LevelDataFAB as LDFAB
from LevelDataFB import LevelDataFluxBox as LDFB
from BoxLayout import BoxLayout as BL
from BoxLayout import DisjointBoxLayout as DBL
from Vects import IntVect,RealVect
from Box import Box
from FAB import FAB
from FluxBox import FluxBox
from FE import Domain 
from functools import wraps
def strip(arg):
    return arg[0]
def id(arg):
    return arg

transmogrifiers = {'bool': strip,
                   'int': strip,
                   'double': strip,
                   'str': strip,
                   'RealVect': RealVect,
                   'IntVect' : IntVect,
                   'id' : id,
                   'FAB': FAB,
                   'FluxBox': FluxBox,
                   'Box': Box,
                   'LevelDataFluxBox': LDFB,
                   'LevelDataFAB': LDFAB,
                   'BoxLayout': BL,
                   'DisjointBoxLayout': DBL, 
                    'Domain' : Domain}

def WhatIs(arg):
    ''' Inspects arg and decide what type it is.
    Tuples should be avoided as native arguments, since in this
    case the transmogrifier will just try to apply the dictionary
    correspoding to the last element of the tuple to the tuple itself.'''


    if isinstance(arg, (str, FAB, FluxBox, Box, LDFAB, LDFB, int, bool, float, list, BL, DBL, IntVect, RealVect, Domain)):
        return 'id'

    try:
        
        return arg[-1]

    except:
        return 'id'




def SOMAR(func):
    @wraps(func)
    def wrapper(*args, **kargs):
        try:
            args = tuple([transmogrifiers[WhatIs(a)](a) for a in args])
        except:
            print([WhatIs(a) for a in args])
            raise(ValueError(""))
        return func(*args, **kargs)
    return wrapper


