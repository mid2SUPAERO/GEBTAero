# coding=UTF-8
import numpy as np
from gebtaero import *
import math

test1 = True
# ~ test1 = False

test2 = True
# ~ test2 = False

test3 = True
# ~ test3 = False

test4 = True
# ~ test4 = False   

test5 = True
# ~ test5= False   

test6 = True
# ~ test6= False   


verbosity = 0
# ~ verbosity = 1

Chord = 1.8288
parameterA = -0.34
EIg2 = 9.77e6
#~ EIg3 = 4e6
EIg3 = 0
GJ = 0.99e6
I11 = 8.64
Nu = 0.18288
Mu = 35.71
SectionLength = 6.096
ndiv = 10

RefPoint = np.zeros([3])
WingAxis = np.zeros([3])
WingAxis[1] = 1.

#~ WingTwist = 0.25*np.pi
WingTwist = 0.

#~ Use of Cross SectionClass
CS = CrossSection()
CS.SetFlexibilityMatrixByIsotropicValues(EIg2,EIg3,GJ)
CS.SetMassMatrix(Mu,1e-4*I11,I11,0.,Nu)

#~ Use of Frame Class
Frame1 = Frame(WingAxis,WingTwist)

#~ Use of WingSection class
GolandSection = WingSection(Chord,parameterA,ndiv,SectionLength,CS,Frame1)

#~ Use of Wing class 
Goland = Wing("Goland",RefPoint)
Goland.AppendWingSection(GolandSection)

#~ Use of Simu class
Vmin = 10
Vmax = 200
Vstep = 5
DeltaV= 0.01
Rho = 1.225
AlphaAC = 0
BetaAC = 0
Ksitol = 1e-6
Lifttol = 1e-3
ModesToPlot = 5
ModesToCompute = 40
Nev = 20
AeroFlag = 3
FreqLim=100

Simu = Simulation(Goland)


if test1:
    #test 1 : undeformed wing flutter speed
    print("test 1 : undeformed wing flutter speed")
    AeroFlag = 1
    Vflutter = Simu.ModalFlutterSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToPlot,ModesToPlot,AlphaAC,BetaAC,Ksitol,verbosity=verbosity)
    assert(math.isclose(Vflutter[0],35.26,abs_tol=1e-1)),"Error flutter speed AeroFlag = 1, flutter speed ="+str(Vflutter[0])+" instead of "+str(35.26)
    Vflutter = Simu.ModalFlutterSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,FreqLim,AlphaAC,BetaAC,KsiObj=1e-6,verbosity=verbosity,arpack=0)
    assert(math.isclose(Vflutter[0],35.26,abs_tol=1e-1)),"Error flutter speed AeroFlag = 1, flutter speed ="+str(Vflutter[0])+" instead of "+str(35.26)
    Vflutter = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=verbosity,mode=1)
    assert(math.isclose(Vflutter[0],35.2,abs_tol=1e-1)),"Error flutter speed AeroFlag = 1, flutter speed ="+str(Vflutter[0])+" instead of "+str(35.26)


    AeroFlag = 2
    Vflutter = Simu.ModalFlutterSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToPlot,ModesToPlot,AlphaAC,BetaAC,Ksitol,verbosity=verbosity)
    assert(math.isclose(Vflutter[0],64.34,abs_tol=1e-1)),"Error flutter speed AeroFlag = 2"
    Vflutter = Simu.ModalFlutterSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,FreqLim,AlphaAC,BetaAC,KsiObj=1e-6,verbosity=verbosity)
    assert(math.isclose(Vflutter[0],64.34,abs_tol=1e-1)),"Error flutter speed AeroFlag = 2"
    Vflutter = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=verbosity,mode=1)
    assert(math.isclose(Vflutter[0],64.34,abs_tol=1e-1)),"Error flutter speed AeroFlag = 2"

    AeroFlag = 3
    Vflutter = Simu.ModalFlutterSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToPlot,ModesToPlot,AlphaAC,BetaAC,Ksitol,verbosity=verbosity)
    assert(math.isclose(Vflutter[0],136.5,abs_tol=1e-1)),"Error flutter speed AeroFlag = 3"
    Vflutter = Simu.ModalFlutterSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,FreqLim,AlphaAC,BetaAC,KsiObj=1e-6,verbosity=verbosity)
    assert(math.isclose(Vflutter[0],136.5,abs_tol=1e-1)),"Error flutter speed AeroFlag = 3"
    Vflutter = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=verbosity,mode=1)
    assert(math.isclose(Vflutter[0],136.5,abs_tol=1e-1)),"Error flutter speed AeroFlag = 3"
    print("#####OK#####")
  
if test2:
    #test 2 : temporal divergence speed
    print("test 2 : temporal divergence speed")  
    DeltaV = 0.01
    Vmin = 5
    Vmax = 300
    VStep = 1.
    Rho = 1.225
    AeroFlag = 1
    DivSpeed = Simu.TemporalDivergenceSpeed(Rho,Vmin,Vmax,VStep,DeltaV,AeroFlag,BetaAC,verbosity=verbosity)
    assert(math.isclose(DivSpeed,252.5,abs_tol=1)),"Error temporal divergence speed"
    AeroFlag = 2
    DivSpeed = Simu.TemporalDivergenceSpeed(Rho,Vmin,Vmax,VStep,DeltaV,AeroFlag,BetaAC,verbosity=verbosity)
    assert(math.isclose(DivSpeed,252.5,abs_tol=1)),"Error temporal divergence speed"
    AeroFlag = 3
    DivSpeed = Simu.TemporalDivergenceSpeed(Rho,Vmin,Vmax,VStep,DeltaV,AeroFlag,BetaAC,verbosity=verbosity)
    assert(math.isclose(DivSpeed,252.5,abs_tol=1)),"Error temporal divergence speed"
    print("#####OK#####")
    
if test3:
    #test 3 : modal divergence speed
    print("test 3 : modal divergence speed")  
    DeltaV = 0.1
    Vmax = 300
    Rho = 1.225
    ModesToPlot = 5
    AeroFlag = 2
    DivSpeed = Simu.ModalDivergenceSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToPlot,ModesToPlot,AlphaAC,BetaAC,CorrCoef=0.9,GravFlag=0,verbosity=verbosity)[0]
    assert(math.isclose(DivSpeed,253,abs_tol=1)),"Error modal divergence speed"
    DivSpeed = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=verbosity,mode=2)
    assert(math.isclose(DivSpeed[0],253,abs_tol=1)),"Error modal divergence speed"
    print("#####OK#####")
    
if test4:
    #test 4 : sweep angle effect
    print("test 4 : sweep angle effect")    
    AeroFlag=3
    #sweep angle:
    Gamma = 15
    Gamma = Gamma*np.pi/180.
    #dihedral angle
    Dihedral = 0.
    Dihedral = Dihedral*np.pi/180.
    # sideslip balance the sweep angle
    BetaAC = Gamma
    WingAxis2 = np.zeros([3])
    WingAxis2[0] = -math.sin(Gamma)
    WingAxis2[1] = math.cos(Gamma)*math.cos(Dihedral)
    WingAxis2[2] = -math.sin(Dihedral)
    WingTwist = 0.

    Frame2 = Frame(WingAxis2,WingTwist)
    GolandSectionSW = WingSection(Chord,parameterA,ndiv,SectionLength,CS,Frame2)
    GolandSW = Wing("GolandSW",RefPoint)
    GolandSW.AppendWingSection(GolandSectionSW)
    SimuSW = Simulation(GolandSW)
    # ~ VecNul=np.zeros(3)
    # ~ SimuSW.Input = InputFile("truc",3,AeroFlag,0,1000,100,0,4,VecNul,0,VecNul,0,SimuSW.Wing,10,Rho,AlphaAC,BetaAC,0.,0.)       
    # ~ SimuSW.Input.WriteInputFile()
    Vflutter = SimuSW.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=1,verbosity=verbosity,mode=1)
    assert(math.isclose(Vflutter[0],136.5,abs_tol=1e-1)),"Error flutter speed swept wing"
    print("#####OK#####")

if test5:
    #test 5 : dihedral angle effect
    print("test 5 : dihedral angle effect")    
    AeroFlag=3
    #sweep angle:
    Gamma = 0.
    Gamma = Gamma*np.pi/180.
    #dihedral angle
    Dihedral = 15
    Dihedral = Dihedral*np.pi/180.
    # sideslip balance the sweep angle
    BetaAC = Gamma
    WingAxis3 = np.zeros([3])
    WingAxis3[0] = -math.sin(Gamma)
    WingAxis3[1] = math.cos(Gamma)*math.cos(Dihedral)
    WingAxis3[2] = -math.sin(Dihedral)
    WingTwist = 0.
    Frame3 = Frame(WingAxis3,WingTwist)
    GolandSectionDA = WingSection(Chord,parameterA,ndiv,SectionLength,CS,Frame3)
    GolandDA = Wing("GolandDA",RefPoint)
    GolandDA.AppendWingSection(GolandSectionDA)
    SimuDA = Simulation(GolandDA)
    Vflutter = SimuDA.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=1,verbosity=verbosity,mode=1)
    assert(math.isclose(Vflutter[0],136.5,abs_tol=1e-1)),"Error flutter speed swept wing"
    print("#####OK#####")
    
if test6:
    #test 6 : twist angle effect
    print("test 6 : twist angle effect")    
    AeroFlag=3
    #sweep angle:
    Gamma = 0.
    Gamma = Gamma*np.pi/180.
    #dihedral angle
    Dihedral = 0.
    Dihedral = Dihedral*np.pi/180.
    # sideslip balance the sweep angle
    BetaAC = Gamma
    WingAxis4 = np.zeros([3])
    WingAxis4[0] = -math.sin(Gamma)
    WingAxis4[1] = math.cos(Gamma)*math.cos(Dihedral)
    WingAxis4[2] = -math.sin(Dihedral)
    WingTwist = np.pi/4.
    AlphaAC = -WingTwist
    Frame4 = Frame(WingAxis4,WingTwist)
    GolandSectionWT = WingSection(Chord,parameterA,ndiv,SectionLength,CS,Frame4)
    GolandWT = Wing("GolandWT",RefPoint)
    GolandWT.AppendWingSection(GolandSectionWT)
    SimuWT = Simulation(GolandWT)
    Vflutter = SimuWT.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=verbosity,mode=1)
    assert(math.isclose(Vflutter[0],136.5,abs_tol=1e-1)),"Error flutter speed swept wing"
    print("#####OK#####")
