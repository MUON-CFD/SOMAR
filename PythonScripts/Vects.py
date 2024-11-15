


class Vect:
    def __init__(self):
        self.Type=None
        self.Vect=None
    
    def __str__(self):
        string = self.Type + " " + str(self.Vect)
        return string

    def toCPP(self):
        print(tuple(i for i in self.Vect))
        return tuple(i for i in self.Vect) 

    def __iter__(self):
        return iter(self.Vect)
        

    def __next__(self):
        return next(self.Vect)

    def __getitem__(self,key):
        return self.Vect[key]

    def __setitem__(self,key,value):
        self.Vect[key]=value
    
    def __len__(self):
        return len(self.Vect)

class IntVect(Vect):
    def __init__(self, Arg):
        if Arg[-1]!="IntVect":
            raise ValueError("IntVect called with non IntVect args")
        self.Vect=tuple(int(i) for i in Arg[0:-1])
        self.Type='IntVect'
    
    

class RealVect(Vect):
    def __init__(self, Arg):
        if Arg[-1]!="RealVect":
            raise ValueError("RealVect called with non RealVect args")
        self.Vect=tuple(float(i) for i in Arg[0:-1])
        self.Type='RealVect'
    
    