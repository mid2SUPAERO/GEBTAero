# coding=UTF-8
import numpy as np
from gebtaero import *
# ~ from time import time

# ~ Start = time()

Chord = 0.05
Ep = 0.001
OffsetY = 0.
OffsetZ = 0.

E = 73e9
Nu = 0.3
G = E/(2+2*Nu)
Rho = 2.7e3
#analytical value for the flexibility matrix 
S = Chord*Ep
#Roark formula p 406 (Chord >>Ep)
a = 0.5*Chord
b = 0.5*Ep
J = a*b**3*(16./3-3.36*b/a*(1-b**4/(12*a**4)))
Ig2 = (Chord*Ep**3)/12
Ig3 = (Chord**3*Ep)/12
S11 = (E*S)**(-1)
S44 = (G*J)**(-1) 
S55 = (E*Ig2)**(-1)
S66 = (E*Ig3)**(-1)
S15 = -S55*OffsetZ
S16 = S66*OffsetY
S11 = S11 + ((OffsetY*S16)**2)**(0.5) + ((OffsetZ*S15)**2)**(0.5)
#analytical value for the mass matrix 
Mu = Rho*S
I22 = Rho*Ig2+Mu*OffsetZ**2
I33 = Rho*Ig3+Mu*OffsetY**2
I11 = I22+I33
I23 = -Mu*OffsetY*OffsetZ

Alu = IsoMaterial(E,Nu,Rho)
Plate = CompositePlate(Chord,OffsetY,OffsetZ)
Plate.AppendPly(CompositePly(Alu,Ep,0))
# ~ Plate.AppendPly(CompositePly(Alu,0.25*Ep,0))
# ~ Plate.AppendPly(CompositePly(Alu,0.25*Ep,0))
# ~ Plate.DisplaySectionDeformation("he20r",1,10,5,2)
CS = CrossSection()
CS.SetFlexibilityMatrixByPlate(Plate,"he20r",1,10,10,RigidX=False,RigidZ=False)
CS.SetMassMatrixByPlate(Plate)
RemoveMeshFiles()
Flex = CS.GetFlexibilityMatrix()
Mass = CS.GetMassMatrix()

# ~ print(Flex)
# ~ print(Mass)

test1 = True
# ~ test1 = False

test2 = True
# ~ test2= False

# ~ End = time()    
# ~ print("Elapsed time : {0}".format(End-Start))

if test1:
    #test 1 : alu plate flexibility matrix versus analytic
    print("test 1 : alu plate flexibility matrix versus analytic")  
    assert(math.isclose(Flex[0,0],S11,rel_tol=5e-2)),"Error 1/ES"
    assert(math.isclose(Flex[3,3],S44,rel_tol=5e-2)),"Error 1/GJ"
    assert(math.isclose(Flex[4,4],S55,rel_tol=1e-2)),"Error 1/EIg2"
    assert(math.isclose(Flex[5,5],S66,rel_tol=1e-2)),"Error 1/EIg3"
    assert(math.isclose(Flex[0,3],0.,abs_tol=1e-6)),"Error non diag coef"
    assert(math.isclose(Flex[0,4],S15,rel_tol=1e-2)),"Error non diag coef"
    assert(math.isclose(Flex[0,5],S16,rel_tol=1e-2)),"Error non diag coef"
    assert(math.isclose(Flex[3,4],0.,abs_tol=1e-6)),"Error non diag coef"
    assert(math.isclose(Flex[3,5],0.,abs_tol=1e-6)),"Error non diag coef"
    assert(math.isclose(Flex[4,5],0.,abs_tol=1e-6)),"Error non diag coef"
    print("#####OK#####")
    
if test2:    
    #test 1 : alu plate mass matrix versus analytic
    print("test 1 : alu plate mass matrix versus analytic")  
    assert(math.isclose(Mass[0,0],Mu,rel_tol=1e-2)),"Error in mu"
    assert(math.isclose(Mass[1,1],Mu,rel_tol=1e-2)),"Error in mu"
    assert(math.isclose(Mass[2,2],Mu,rel_tol=1e-2)),"Error in mu"
    assert(math.isclose(Mass[3,3],I11,rel_tol=1e-2)),"Error in i1"
    assert(math.isclose(Mass[4,4],I22,rel_tol=1e-2)),"Error in i2"
    assert(math.isclose(Mass[5,5],I33,rel_tol=1e-2)),"Error in i3"
    assert(math.isclose(Mass[0,3],0.,abs_tol=1e-12)),"Error non diag coef"
    assert(math.isclose(Mass[0,4],-Mu*OffsetZ,rel_tol=1e-2)),"Error non diag coef"
    assert(math.isclose(Mass[0,5],Mu*OffsetY,rel_tol=1e-2)),"Error non diag coef"
    assert(math.isclose(Mass[1,3],Mu*OffsetZ,rel_tol=1e-2)),"Error non diag coef"
    assert(math.isclose(Mass[2,3],-Mu*OffsetY,rel_tol=1e-2)),"Error non diag coef"
    assert(math.isclose(Mass[3,4],0.,abs_tol=1e-12)),"Error non diag coef"
    assert(math.isclose(Mass[3,5],0.,abs_tol=1e-12)),"Error non diag coef"
    assert(math.isclose(Mass[4,5],I23,rel_tol=1e-2)),"Error non diag coef"
    print("#####OK#####")
