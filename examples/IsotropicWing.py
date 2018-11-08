## This example show a number of applications on isotropic wing using data from literature test cases (the example uses the "Patil wing")
#after the problem setting each pplication is in a if codeblock. To use an application, decomment "if True:" and comment "if False:"

import numpy as np
import os
from gebtaero import *

# Wing parameters
Chord = 1.      #Wing chord
parameterA = 0      #Distance between mid chord and reference line (cf article)
EIg2 = 2e4      #Bending stiffness (span-wise)
EIg3 = 0        #Bending stiffness (chord-wise) ( if 0 : the ddl is suppressed)
GJ = 1e4        #Torsional stiffness
I22 = 1e-5  #Mass moment of inertia around chord axis (<<I33)
I33 = 0.1-I22      #Mass moment of inertia around reference line
I23 = 0.        #product of inertia
Nu = 0      # Distance between  Center of Gravity and Elastic Axis (>0 if the CG is behind the EA)
Mu = 0.75   #Mass per unit length
SectionLength = 16  #Wing section length (Wing length in case of a unique section)
ndiv = 24   #Number of finite element per wing section
RefPoint = np.zeros([3])        #Wing origin
WingAxis = np.zeros([3])        #Wing axis
WingAxis[1] = 1.        # the wing axis is ya (no sweep angle and no dihedral)
WingTwist = 0.      # no twist angle

#~ Cross Section definition
CS = CrossSection()
CS.SetFlexibilityMatrixByIsotropicValues(EIg2,EIg3,GJ)
CS.SetMassMatrix(Mu,I22,I33,I23,Nu)

# Frame definition (frame b)
Frame = Frame(WingAxis,WingTwist)

# Creation of a wing section
Section = WingSection(Chord,parameterA,ndiv,SectionLength,CS,Frame)

# Creation of a wing
IsoWing = Wing("Patil",RefPoint)
IsoWing.AppendWingSection(Section)

# Creation of a simulation object
Simu = Simulation(IsoWing)

# algorithm parameters
verbosity = 0     # 0 = without log, 1 = with log

#~ Flight parameters
AlphaAC = 0     # Aircraft angle of attack
BetaAC = 0      # Aircraft yaw angle
Rho = 0.0889     # Air density
Vinf = 10        # Upstream Velocity
AeroFlag = 3        # type of aerodynamic model used : 0 = no aero, 1 = quasi-steady, 2= quasi-steady with added mass, 3 = unsteady (Peters)
GravFlag = 1        # 0 = without gravity forces, 1 = with gravity forces

###########################
#---------------------------------------------------------#
# static wing structural modes (no velocity)
# ~ if True:
if False:    
    Nev = 8      #maximum number of modes to be computed 
    Vinf = 0.       # Upstream Velocity
    vtk = 1     # 0 = no output file, 1 = creation of a folder with output vtk file (vtk format to be used in paraview)
    Modes = Simu.Eigenvalues(Vinf,Rho,AlphaAC,BetaAC,AeroFlag,Nev,GravFlag=0,verbosity=1,vtk=1)
    VtkFolder = os.getcwd()+"/PatilModes.dat_vtk/"
    VtkName = "PatilModes.dat"
    GebtPlot.WriteParaviewScript(VtkFolder,VtkName)
    #Launch the paraview script
    RunParaviewScript("pvscript.py")
###########################
#---------------------------------------------------------#
# static loads due to gravity at the wing root
# ~ if True:
if False:   
    Static = Simu.StaticLoads(Vinf,Rho,AlphaAC,BetaAC,GravFlag=1,verbosity=1)

###########################
#---------------------------------------------------------#
# Computation of critical speed (divergence and flutter)
# ~ if True:
if False:   
    DeltaV = 0.1        # Critical velocity tolerance
    Vstep = 1       # Velocity step used in the algorithm
    mode = 3        # Algorithm setting : 0 = the first critical speed, 1 = flutter only, 2 = divergence only, 3 =  flutter and divergence
    Vmin = 0.       # Lower boundary of the search interval
    Vmax = 100.       # Upper boundary of the search interval
    Vcritique = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=GravFlag,verbosity=verbosity,mode=mode)
    print("Flutter speed =",Vcritique[0])
    print("Flutter frequency =",Vcritique[1])
    print("Divergence speed =",Vcritique[2])

###########################
#---------------------------------------------------------#
# Plotting of the first aeroelastic modes
# ~ if True:
if False:   
    Nstep = 50    #Number of frequency values in the interval
    ModesToPlot = 5     #Number of Aeroelastic modes to plot
    Vmin = 0.       # Lower boundary of the plotting interval
    Vmax = 50.       # Upper boundary of the plotting interval
    DeltaV = 0.01       # Velocity tolerance of the mode correlation algorithm
    CorrCoef=0.99       # Minimal mode correlation coefficient allowed
    Modes = Simu.EigenTabSorted(Rho,Vmax,DeltaV,Nstep,AeroFlag,ModesToPlot,ModesToPlot,AlphaAC,BetaAC,CorrCoef=0.99,GravFlag=GravFlag,verbosity=1)
    #ReducedDamping: True if the real part of the mode is converted into reduced damping; DampAxis: True if the damping scale is reduced to Damp Scale value
    GebtPlot.EigenFreqDamping(Modes[0],Modes[1],ReducedDamping=True,DampAxis=True,DampScale=0.01)      


###########################
#---------------------------------------------------------#
# making a temporal simulation of the flutter phenomena
if True:
# ~ if False:   
    # a modal flutter speed and frequency is first calculated, 
    Vmin = 0.
    Vmax = 100
    DeltaV = 0.1    
    NbPeriod = 100      #number of period used to calculate the initial simulation time; the simulation is adjusted after a first simulation with a shorter period to catch the phenomena before too much displacement
    StepByPeriod = 100      #number of time step per period
    CoefPerturb = 0.0001        #a vertical speed to initiate the flutter is applied during one period with a coefficient relative to the upstream velocity
    CoefVinf = 1.2    #the simulation is done a V = ModalFlutterVelocity*CoefVinf
    Nvtk = 1000     #number of vtk output file
    FlutterLimit = 0.4      #value of the maximal angular deformation (radian) of the wing which will trigger the temporal flutter state
    # Realisation of temporal simulation and creation a a vtk file folder "PatilTempFlutter.dat_vtk/" in the work directory 
    Simu.FlutterVtk(Rho,Vmin,Vmax,DeltaV,AeroFlag,AlphaAC,BetaAC,NbPeriod,StepByPeriod,CoefPerturb,CoefVinf,Nvtk,FlutterLimit,GravFlag=GravFlag,verbosity=1)
    # write a paraview script with the vtk folder
    VtkFolder = os.getcwd()+"/PatilTempFlutter.dat_vtk/"
    VtkName = "PatilTempFlutter.dat"
    GebtPlot.WriteParaviewScript(VtkFolder,VtkName)
    #Launch the paraview script
    RunParaviewScript("pvscript.py")



