# coding=UTF-8
import numpy as np
from gebtaero import *
import math

test1 = True
# ~ test1 = False

test2 = True
# ~ test2 = False

# ~ verbosity = 0
verbosity = 1

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
Vmin = 1
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
FreqLim=100.

Simu = Simulation(Goland)

if test1:
    #test 1 : temporal goland flutter speed
    print("test 1 : temporal goland flutter speed")
    NbPeriod = 1000
    StepByPeriod = 10
    CoefPerturb = 1e-6
    DeltaV = 1
    Vmin = 30
    VFlutterTemp = Simu.TemporalFlutterSpeed(Rho,Vmin,Vmax,DeltaV,AeroFlag,AlphaAC,BetaAC,NbPeriod,StepByPeriod,CoefPerturb,GravFlag=0,FlutterLimit=0.1,verbosity=verbosity)
    assert(math.isclose(VFlutterTemp,137.,rel_tol=0.02)),"Error temporal flutter speed"
    print("#####OK#####")
    
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


#~ Use of Cross SectionClass
CS2 = CrossSection()
CS2.SetFlexibilityMatrixByIsotropicValues(EIg2,EIg3,GJ)
CS2.SetMassMatrix(Mu,1e-3*I11,I11,0.,Nu)

#~ Use of WingSection class
PatilSection = WingSection(Chord,parameterA,ndiv,SectionLength,CS2,Frame1)

#~ Use of Wing class 
Patil = Wing("Patil2",RefPoint)
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

if test2:
    #test 2 : temporal patil flutter speed
    print("test 2 : temporal patil flutter speed")
    NbPeriod = 1000
    StepByPeriod = 10
    CoefPerturb = 1e-6
    DeltaV = 0.5
    Nev = 60
    VFlutterTemp = Simu.TemporalFlutterSpeed(Rho,Vmin,Vmax,DeltaV,AeroFlag,AlphaAC,BetaAC,NbPeriod,StepByPeriod,CoefPerturb,FlutterLimit=0.1,verbosity=verbosity)
    assert(math.isclose(VFlutterTemp,32.,abs_tol=1)),"Error temporal flutter speed"
    print("#####OK#####")
    
    
