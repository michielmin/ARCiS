import multiprocessing 

import matplotlib.pyplot as plt
import numpy as np
import pyARCiS

class A:
    def __init__(self):
        pyARCiS.pyinit("inputSimple.dat","output")
        self.lam=np.empty([pyARCiS.pyex.nlam],dtype=float)
        self.trans=np.empty([pyARCiS.pyex.nlam],dtype=float)
        self.emis=np.empty([pyARCiS.pyex.nlam],dtype=float)
        self.P=np.empty([pyARCiS.pyex.nr],dtype=float)
        self.T=np.empty([pyARCiS.pyex.nr],dtype=float)
    def setpar(self, Rp):
        pyARCiS.pysetvalue('Rp',Rp)
        self.Rp=Rp
    def getspec(self,i):
        pyARCiS.pycomputemodel()
        pyARCiS.pygettrans(self.lam,self.trans)
        return self.Rp,self.trans

def pool_init(i):
    global a
    a = A()

def pool_getspec(*args):
    global a
    Rp, trans = a.getspec(*args)
    return Rp, trans
    
def pool_setpar(*args):
    global a
    a.setpar(*args)
    
if __name__ == "__main__":
    #a = A() #no point in instantiating A in main process as the children will have separate instances
    p1=multiprocessing.Pool(4)

    p1.map(pool_init,[1,2,3,4])
    p1.map(pool_setpar,[0.15,0.25,0.35,0.45])
    res=p1.map(pool_getspec,[1,2,3,4])
    for r in res:
    	Rp,trans=r
    	print(Rp,trans)


