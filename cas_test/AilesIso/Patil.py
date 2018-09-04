# coding=UTF-8
import numpy as np
from gebtaero import *
import math

# ~ test1 = True
test1 = False

# ~ test2 = True
test2 = False

# ~ test3 = True
test3 = False

test4 = True
# ~ test4 = False

# ~ test5 = True
test5 = False

verbosity = 0
# ~ verbosity = 1

Chord = 1.
parameterA = 0
EIg2 = 2e4
EIg3 = 0
GJ = 1e4
I11 = 0.1
Nu = 0.
Mu = 0.75
SectionLength = 16.
ndiv = 10

RefPoint = np.zeros([3])
WingAxis = np.zeros([3])
WingAxis[1] = 1.

#~ WingTwist = 0.25*np.pi
WingTwist = 0.

#~ Use of Cross SectionClass
CS = CrossSection()
CS.SetFlexibilityMatrixByIsotropicValues(EIg2,EIg3,GJ)
CS.SetMassMatrix(Mu,1e-3*I11,I11,0.,Nu)

#~ Use of Frame Class
Frame1 = Frame(WingAxis,WingTwist)

#~ Use of WingSection class
PatilSection = WingSection(Chord,parameterA,ndiv,SectionLength,CS,Frame1)

#~ Use of Wing class 
Patil = Wing("Patil",RefPoint)
Patil.AppendWingSection(PatilSection)

#~ Use of Simu class
Vmin = 0
Vmax = 50
Vstep = 5
DeltaV= 0.1
Rho = 0.0889
AlphaAC = 0
BetaAC = 0
Ksitol = 1e-1
Ksiobj = 1e-6
Lifttol = 1e-3
ModesToPlot = 4
ModesToCompute = 120
Nev = 20
AeroFlag = 3
FreqLim = 50

Simu = Simulation(Patil)

if test1:
    #test 1 : undeformed wing flutter speed
    print("test1 : undeformed wing flutter speed")
    AeroFlag = 1
    Nu = -0.1
    ModesToPlot = 4
    Nev = 20
    CS.SetMassMatrix(Mu,1e-4*I11,I11,0.,Nu)
    Vflutter = Simu.ModalFlutterSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,FreqLim,AlphaAC,BetaAC,KsiObj=1e-6,verbosity=verbosity)
    assert(math.isclose(Vflutter[0],14.83,abs_tol=1e-1)),"Error flutter speed Nu = -0.1, AeroFlag = 1"
    Vflutter = Simu.ModalFlutterSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToPlot,ModesToPlot,AlphaAC,BetaAC,Ksiobj,verbosity=verbosity)
    assert(math.isclose(Vflutter[0],14.83,abs_tol=1e-1)),"Error flutter speed Nu = -0.1, AeroFlag = 1"
    Vflutter = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=verbosity,mode=1)
    assert(math.isclose(Vflutter[0],14.83,abs_tol=1e-1)),"Error flutter speed Nu = -0.1, AeroFlag = 1"

    AeroFlag = 2
    Nu = -0.1
    CS.SetMassMatrix(Mu,1e-3*I11,I11,0.,Nu)
    Vflutter = Simu.ModalFlutterSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,FreqLim,AlphaAC,BetaAC,KsiObj=1e-6,verbosity=verbosity)
    assert(math.isclose(Vflutter[0],25.85,abs_tol=1e-1)),"Error flutter speed Nu = -0.1, AeroFlag = 2"
    Vflutter = Simu.ModalFlutterSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToPlot,ModesToPlot,AlphaAC,BetaAC,Ksiobj  ,verbosity=verbosity)
    assert(math.isclose(Vflutter[0],25.85,abs_tol=1e-1)),"Error flutter speed Nu = -0.1, AeroFlag = 2"
    Vflutter = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=verbosity,mode=1)
    assert(math.isclose(Vflutter[0],25.85,abs_tol=1e-1)),"Error flutter speed Nu = -0.1, AeroFlag = 2"

    AeroFlag = 1
    Nu = -0.05
    CS.SetMassMatrix(Mu,1e-3*I11,I11,0.,Nu)
    Vflutter = Simu.ModalFlutterSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,FreqLim,AlphaAC,BetaAC,KsiObj=1e-6,verbosity=verbosity)
    assert(math.isclose(Vflutter[0],6.32,abs_tol=1e-1)),"Error flutter speed Nu = -0.05, AeroFlag = 1"
    Vflutter = Simu.ModalFlutterSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToPlot,ModesToPlot,AlphaAC,BetaAC,Ksiobj,verbosity=verbosity)
    assert(math.isclose(Vflutter[0],6.32,abs_tol=1e-1)),"Error flutter speed Nu = -0.05, AeroFlag = 1"
    Vflutter = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=verbosity,mode=1)
    assert(math.isclose(Vflutter[0],6.16,abs_tol=1e-1)),"Error flutter speed Nu = -0.05, AeroFlag = 1"

    AeroFlag = 2
    Nu = -0.05
    CS.SetMassMatrix(Mu,1e-3*I11,I11,0.,Nu)
    Vflutter = Simu.ModalFlutterSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,FreqLim,AlphaAC,BetaAC,KsiObj=1e-6,verbosity=verbosity)
    assert(math.isclose(Vflutter[0],14.83,abs_tol=2e-1)),"Error flutter speed Nu = -0.05, AeroFlag = 2"
    Vflutter = Simu.ModalFlutterSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToPlot,ModesToPlot,AlphaAC,BetaAC,Ksiobj,verbosity=verbosity)
    assert(math.isclose(Vflutter[0],14.83,abs_tol=1e-1)),"Error flutter speed Nu = -0.05, AeroFlag = 2"
    Vflutter = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=verbosity,mode=1)
    assert(math.isclose(Vflutter[0],14.83,abs_tol=1e-1)),"Error flutter speed Nu = -0.05, AeroFlag = 2"

    AeroFlag = 3
    Nu = 0.
    Nev = 40
    CS.SetMassMatrix(Mu,1e-3*I11,I11,0.,Nu)
    Vflutter = Simu.ModalFlutterSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,FreqLim,AlphaAC,BetaAC,KsiObj=1e-6,verbosity=verbosity)
    assert(math.isclose(Vflutter[0],32.12,abs_tol=1e-1)),"Error flutter speed Nu = 0, AeroFlag = 3"
    Vflutter = Simu.ModalFlutterSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToPlot,ModesToPlot,AlphaAC,BetaAC,Ksiobj,verbosity=verbosity)
    assert(math.isclose(Vflutter[0],32.12,abs_tol=1e-1)),"Error flutter speed Nu = 0, AeroFlag = 3"
    Vflutter = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=verbosity,mode=1)
    assert(math.isclose(Vflutter[0],32.12,abs_tol=1e-1)),"Error flutter speed Nu = 0, AeroFlag = 3"
    print("#####OK#####")

if test2:
    Nev = 60
    ModesToPlot = 6
    FreqLim = 20

    #test 2 : deformed wing flutter speed
    print("test 2 : deformed wing flutter speed")
    VflutterDef = Simu.ModalFlutterSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,FreqLim,AlphaAC,BetaAC,KsiObj=1e-6,GravFlag=1,verbosity=verbosity)
    assert(math.isclose(VflutterDef[0],23.35,abs_tol=2e-1)),"Error Vflutter deformed"
    VflutterDef = Simu.ModalFlutterSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToPlot,ModesToPlot,AlphaAC,BetaAC,Ksiobj,GravFlag=1,verbosity=verbosity)
    assert(math.isclose(VflutterDef[0],23.35,abs_tol=2e-1)),"Error Vflutter deformed"
    VflutterDef = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=1,verbosity=verbosity,mode=1)
    assert(math.isclose(VflutterDef[0],23.35,abs_tol=2e-1)),"Error Vflutter deformed"
    print("#####OK#####")

if test3:
    #test 3 :  divergence speed
    print("test 3 :  divergence speed")  
    DeltaV = 0.01
    Vmin = 5
    Vmax = 100
    Vstep = 1.
    AeroFlag = 1
    Rho = 0.0889
    DivSpeed = Simu.TemporalDivergenceSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,BetaAC,verbosity=verbosity)
    assert(math.isclose(DivSpeed,37.25,abs_tol=1e-1)),"Error temporal flutter speed"
    DivSpeed = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=verbosity,mode=2)[0]
    assert(math.isclose(DivSpeed,37.25,abs_tol=1e-1)),"Error temporal flutter speed"
    Rho = 1.225
    DivSpeed = Simu.TemporalDivergenceSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,BetaAC,verbosity=verbosity)
    assert(math.isclose(DivSpeed,10.01,abs_tol=1e-1)),"Error temporal flutter speed"
    DivSpeed = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=verbosity,mode=2)[0]    
    assert(math.isclose(DivSpeed,10.01,abs_tol=1e-1)),"Error temporal flutter speed"
    AeroFlag = 2
    Rho = 0.0889
    DivSpeed = Simu.TemporalDivergenceSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,BetaAC,verbosity=verbosity)
    assert(math.isclose(DivSpeed,37.25,abs_tol=1e-1)),"Error temporal flutter speed"
    DivSpeed = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=verbosity,mode=2)[0]    
    assert(math.isclose(DivSpeed,37.25,abs_tol=1e-1)),"Error temporal flutter speed"
    Rho = 1.225
    DivSpeed = Simu.TemporalDivergenceSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,BetaAC,verbosity=verbosity)
    assert(math.isclose(DivSpeed,10.01,abs_tol=1e-1)),"Error temporal flutter speed"
    AeroFlag = 3
    Rho = 0.0889
    DivSpeed = Simu.TemporalDivergenceSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,BetaAC,verbosity=verbosity)
    assert(math.isclose(DivSpeed,37.25,abs_tol=1e-1)),"Error temporal flutter speed"
    DivSpeed = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=verbosity,mode=2)[0]    
    assert(math.isclose(DivSpeed,37.25,abs_tol=1e-1)),"Error temporal flutter speed"
    Rho = 1.225
    DivSpeed = Simu.TemporalDivergenceSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,BetaAC,verbosity=verbosity)
    assert(math.isclose(DivSpeed,10.01,abs_tol=1e-1)),"Error temporal flutter speed"
    DivSpeed = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=verbosity,mode=2)[0]    
    assert(math.isclose(DivSpeed,10.01,abs_tol=1e-1)),"Error temporal flutter speed"
    print("#####OK#####")
    
if test4:
    # test 4 : structural modes
    print("test4 : structural modes")
    EIg3 = 4e6
    CS0 = CrossSection()
    CS0.SetFlexibilityMatrixByIsotropicValues(EIg2,EIg3,GJ)
    CS0.SetMassMatrix(Mu,1e-3*I11,I11,0.,Nu)
    ndiv = 20
    PatilSection0 = WingSection(Chord,parameterA,ndiv,SectionLength,CS0,Frame1)
    Patil0 = Wing("Patil0",RefPoint)
    Patil0.AppendWingSection(PatilSection0)
    Simu0 = Simulation(Patil0)
    AeroFlag = 0
    Vinf = 0.
    Nev = 12
    Modes = Simu0.Eigenvalues(Vinf,Rho,AlphaAC,BetaAC,AeroFlag,Nev,verbosity=verbosity)
    Freq=Modes[:,1]
    Freq.sort()
    assert(math.isclose(Freq[0],0.357,abs_tol=1e-2)),"Error mode 1"
    assert(math.isclose(Freq[1],2.257,abs_tol=1e-2)),"Error mode 2"
    assert(math.isclose(Freq[2],4.941,abs_tol=1e-2)),"Error mode 3"
    assert(math.isclose(Freq[3],5.048,abs_tol=1e-2)),"Error mode 4"
    assert(math.isclose(Freq[4],6.428,abs_tol=1e-2)),"Error mode 5"
    print("#####OK#####")
    
if test5:
#test 5 : sweep angle effect
    print("test 5 : sweep angle effect")    
    AeroFlag = 3
    Rho = 0.0889
    ndiv=10
    #sweep angle:
    Gamma = 45
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

    Frame2 = Frame(WingAxis2,WingTwist)
    PatilSectionSW = WingSection(Chord,parameterA,ndiv,SectionLength,CS,Frame2)
    PatilSW = Wing("PatilSW",RefPoint)
    PatilSW.AppendWingSection(PatilSectionSW)
    SimuSW = Simulation(PatilSW)

    Vflutter = SimuSW.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=1,verbosity=verbosity,mode=1)
    assert(math.isclose(Vflutter[0],23.35,abs_tol=2e-1)),"Error sweep angle effect"
    print("#####OK#####")
