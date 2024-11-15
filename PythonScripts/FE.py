
import PyGlue
import numpy as np

class Domain:
    """ defines the FE applied to the whole domain"""
    def __init__(self,Args):
        """ Construct the Domain as passed by SOMAR Args[0] = nodes,
                                                    Args[1] = elements,
                                                    Args[2] = hash"""
        if (Args[-1]!= 'Domain'):
            raise ValueError("Domain called with the wrong arguments")
        
        
        self.nNodes = Args[0][0]
        self.nFaces = Args[1][0]
        self.dim = int(Args[0][1][0])
        self.nVertices = int(Args[1][1][0])

        self.nodes = PyGlue.make_from_memView(Args[0][2], [self.nNodes, self.dim], dtype=np.float64)
        
        self.elements = PyGlue.make_from_memView(Args[1][2], [self.nFaces, self.nVertices], dtype=np.uint64)
       
        self.hash=Args[2][0]
    
    def plot(self):
        """ plot the current domain mesh using pyVista"""
        import pyvista as pv
        
        if self.dim == 3:
            faces = np.zeros((self.nFaces,self.nVertices+1),dtype=np.uint64)
            faces[:,0] = self.nVertices
            faces[:,1:] = self.elements
        
        
            surf = pv.PolyData(self.nodes, faces)
            surf.plot(show_edges=True)
        else:
            print("I only plot 3D stuff")
        
        

        
        