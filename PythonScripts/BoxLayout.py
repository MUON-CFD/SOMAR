from Box import Box as Box

import pickle 
class BoxLayout:

    def __init__(self,Arg):

        self.Boxes, self.PIDs=self._BL(Arg)
    
    def write(self,filename):
        with open(filename, 'wb') as out_file:
            pickle.dump(self, out_file)

    @classmethod
    def read(cls,filename):
        with open(filename, 'rb') as in_file:
            BL=pickle.load(in_file)
        return BL


    def _BL(self,Args):
        """ unpacks a BoxLayout into two tuples containing boxes and ids"""
        
        if Args[-1] != "BoxLayout":
            raise ValueError("BoxLayout called with non BoxLayout Args")

        Boxes=tuple([Box(F) for F in Args[0]])
        PIDs=tuple([p[0] for p in Args[1]])
        
        return (Boxes, PIDs)

    def __str__(self):
        
        string =   "BoxLayout \n"
        string += "=================================\n"
        for Id,box in zip(self.PIDs,self.Boxes):
            string += " Proc Id " + str(Id) + " -- " + str(box) + "\n"
        
        
        return string

    def toCPP(self):
        n=len(self.PIDs)
        output=[n]
        for b,id in zip(self.Boxes, self.PIDs):
            output.append((b.LoEnd()+b.HiEnd()))
            output.append(id)
        
        return tuple(output)

class DisjointBoxLayout(BoxLayout):

    def __init__(self, Arg):
        if Arg[-1] != "DisjointBoxLayout":
            raise ValueError("DisjointBoxLayout called with non DisjointBoxLayout Args")

        super().__init__(Arg[0:2]+('BoxLayout',))
       
        self.Domain=ProblemDomain(Arg[2][0])

    def __str__(self):
        string = "Disjoint"+super().__str__()
        string+= str(self.Domain)
        return string

    def toCPP(self):
        
        return super().toCPP()+(self.Domain.toCPP(),)

    

class ProblemDomain:
    def __init__(self, Arg):
        if Arg[-1] != "ProblemDomain":
            print(Arg[-1])
            raise ValueError("ProblemDomain called with non ProblemDomain Args")

        self.Box=Box(Arg[0])
        self.isPeriodic=Arg[1][0:-1]

    def __str__(self):
        string = "Problem Domain " + str(self.Box)
        return string

    def toCPP(self):
        return self.Box.toCPP()+(tuple(p for p in self.isPeriodic),)


