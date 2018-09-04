# coding=UTF-8
import numpy as np
from gebtaero import *
import math
import os
# ~ from time import time

# ~ Start = time()

E = 73e9
nu = 0.3
G = E/(2+2*nu)
Alu = IsoMaterial(E,0.3,2700.)
Mesh = ExternalMesh(os.path.dirname(os.path.abspath(__file__))+"/mesh_disque/disque.unv",UnvConv=True)
#~ Mesh = ExternalMesh("/mesh_disque/disque.unv",UnvConv=True)
Mesh.AppendComponent("All",Alu)
# create cross section
CS = CrossSection()
CS.SetFlexibilityMatrixByMesh(Mesh,PlaneSection=False,AtElasticCenter=False,RigidX=False,RigidZ=False)
CS.SetMassMatrixByMesh(Mesh,SymY=True,AtElasticCenter=False)
Flex = CS.GetFlexibilityMatrix()
Mass = CS.GetMassMatrix()


test1 = True
#~ test1 = False

test2 = True
#~ test2 = False

# ~ End = time()    
# ~ print("Elapsed time : {0}".format(End-Start))


# ~ S11=1./(E*np.pi)
# ~ S44 = 1./(G*0.5*np.pi)
# ~ S55 = 1./(E*0.25*np.pi)
# ~ S66 = 1./(E*0.25*np.pi)
#relative error %
# ~ print(100*(Flex[0,0]-S11)/S11)
# ~ print(100*(Flex[3,3]-S44)/S44)
# ~ print(100*(Flex[4,4]-S55)/S55)
# ~ print(100*(Flex[5,5]-S66)/S66)

if test1:
    #test 1 : disque alu flexibility matrix versus analytic
    print("test 1 : disque alu flexibility matrix versus analytic")  
    assert(math.isclose(Flex[0,0],4.36e-12,abs_tol=1e-13)),"Error 1/ES"
    assert(math.isclose(Flex[3,3],2.26e-11,abs_tol=1e-12)),"Error 1/GJ"
    assert(math.isclose(Flex[4,4],1.74e-11,abs_tol=1e-12)),"Error 1/EIg2"
    assert(math.isclose(Flex[5,5],1.74e-11,abs_tol=1e-12)),"Error 1/EIg3"
    assert(math.isclose(Flex[0,3],0.,abs_tol=1e-12)),"Error non diag coef"
    assert(math.isclose(Flex[0,4],0.,abs_tol=1e-12)),"Error non diag coef"
    assert(math.isclose(Flex[0,5],0.,abs_tol=1e-12)),"Error non diag coef"
    assert(math.isclose(Flex[3,4],0.,abs_tol=1e-12)),"Error non diag coef"
    assert(math.isclose(Flex[3,5],0.,abs_tol=1e-12)),"Error non diag coef"
    assert(math.isclose(Flex[4,5],0.,abs_tol=1e-12)),"Error non diag coef"
    
    print("#####OK#####")
    
if test2:
    #test 2 : Elastic center position versus analytic
    print("test 2 : Elastic center position versus analytic")  
    EC = ComputeElasticCenterFromFlexMat(Flex)
    assert(math.isclose(EC[0],0,abs_tol=1e-13)),"Error in Elastic center position"
    assert(math.isclose(EC[1],0,abs_tol=1e-13)),"Error in Elastic center position"
    print("#####OK#####")

RemoveFiles()


