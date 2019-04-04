## This example show the construction of a composite plate made of an AIREX foam core surrounded by two carbon/epoxy plies.
# The lower ply has a fix orientation of -45°, the upper ply has a variable orientation.
# critical speed are calculated for each upper ply orientation and plotted
import numpy as np
import matplotlib.pyplot as plt
from gebtaero import *

Chord = 0.2 # width od the composite plate
parameterA = 0  #Distance between mid chord and reference line (cf article)
SectionLength = 0.4 # length of the composite plate
ndiv = 20   # #Number of finite element of the plate

RefPoint = np.zeros([3])    # Plate origin
WingAxis = np.zeros([3])   
WingAxis[1] = 1.    #beam axis

WingTwist = 0.  #beam twisting


# ~ EpComp = 0.25e-3
EpComp = 0.125e-3   #carbon ply thickness
# ~ EpAIREX = 1.2e-3
EpAIREX = 2e-3  #foam core thickness

# ~ PLA = IsoMaterial(2.7e9,0.35,1.25e3)
AIREX = IsoMaterial(66e6,0.32,80) # AIREX C70.75
T700 = OrthoMaterial(148e9,10e9,0.3,4.6e9,1.6e3)    #T700 carbon/epoxy
T300 = OrthoMaterial(181e9,10.3e9,0.28,7.17e9,1.6e3)    #T300 carbon epoxy
FdV = OrthoMaterial(45e9,12e9,0.3,4.5e9,2080)   # glass fiber/epoxy
Alu = IsoMaterial(73e9,0.3,2.7e3)   #Aluminium alloy
CS = CrossSection()     #Cross Section constructor
#~ Use of Frame Class
Frame = Frame(WingAxis,WingTwist)   #Definition of frame b of the plate

#~ Use of WingSection class
AileSection = WingSection(Chord,parameterA,ndiv,SectionLength,CS,Frame) # wing section constructor

#~ Use of Wing class 
Aile = Wing("Aile",RefPoint)    # Wing constructor
Aile.AppendWingSection(AileSection) #append a section to the wing (composite plate)
#~ Use of Simu class
Vmin = 0    # Lower boundary of the search interval
Vmax =100   # Upper boundary of the search interval
Vstep = 1   # Velocity step used in the algorithm
DeltaV= 0.01    # Critical velocity tolerance
Rho = 1.225 # Air density
AlphaAC = 0 # Aircraft angle of attack
BetaAC = 0  # Aircraft yaw angle
AeroFlag = 3    # type of aerodynamic model used : 0 = no aero, 1 = quasi-steady, 2= quasi-steady with added mass, 3 = unsteady (Peters)
Nev = 30    #maximum number of modes to be computed 
Ksitol = 1e-6   # reduced damping tolerance for a critical speed computation
ModesToPlot = 2 # number of modes to follow
Nstep = 50  # number of velocity step to compute for a mode plotting
GravFlag = 1    # 0 = without gravity forces, 1 = with gravity forces

verbosity = 1   # 0 no log output, 1 log output
# ~ verbosity = 0

Simu = Simulation(Aile)
Result = []

# create a list containing the different upper ply orientations
NbOri = 4  # the number of upper ply orientation (equally distributed between -90 and +90)
OriStep = 180/NbOri

#for each orientation, create a plate, compute flexibility and mass matrix to update the cross section parameters and th
for i in range(NbOri+1):
    Ori = OriStep*(i)-90
    Plate = CompositePlate(Chord)   #reset the Plate object
    Plate.AppendPly(CompositePly(T300,EpComp,-45))  #append lower ply
    Plate.AppendPly(CompositePly(AIREX,EpAIREX,0))  #append foam core
    Plate.AppendPly(CompositePly(T300,EpComp,Ori))  #append upper ply
    CS.SetFlexibilityMatrixByPlate(Plate,"he8i",1,10,10,RigidX=True,RigidZ=True)    # compute flexibility matrix (overwrite the previous flexibility matrix). RigidX = True : suppression of traction along beam axis ddl. RigidZ = True : suppression of bending in lag axis ddl
    Flex = CS.GetFlexibilityMatrix()
    CS.SetMassMatrixByPlate(Plate)  #compute flexibility matrix (overwrite the previous flexibility matrix)
    Vflutter = Simu.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=GravFlag,verbosity=verbosity,mode=3) # compute critical speed (both divergence and flutter) for the current orientation
    Result.append([Ori,Vflutter[0],Vflutter[2],1./Flex[3,3],1./Flex[4,4],Vflutter[1]])  # append a result list with upper ply orientation, flutter speed, divergence speed, torsional stiffness, bending stiffness spanwise and flutter frequency

Result = np.array(Result)

# create a plot using the Result numpy array(see matplotlib documentation)
fig, axs = plt.subplots(2,1)
ax1 = axs[0]
ax3 = axs[1]
flutter, = ax1.plot(Result[:,0],Result[:,1],linewidth=3,label='flutter speed',color='green',marker='*',markersize=10,linestyle='solid')
ax1.set_ylabel('flutter speed (m/s)',fontsize=25)
ax1.set_xlim(-90,90)
ax1.tick_params(labelsize=20)

ax2 = ax1.twinx()
freq, = ax2.plot(Result[:,0],Result[:,5],linewidth=3, label='flutter frequency',color='blue',marker='o',markersize=8,linestyle='dashed')
ax2.set_ylabel('frequency (Hz)',fontsize=20)
ax2.tick_params(labelsize=20)
plt.legend(handles=[flutter,freq],loc=7, prop={'size': 25})
plt.grid()

div, = ax3.plot(Result[:,0],Result[:,2],linewidth=3,label='divergence speed',color='orange',marker='*',markersize=10,linestyle='solid')
ax3.set_xlabel(' upper ply orientation (°)',fontsize=25)
ax3.set_ylabel('divergence speed (m/s)',fontsize=25)
ax3.set_xlim(-90,90)
ax3.tick_params(labelsize=20)

ax4 = ax3.twinx()
tors, = ax4.plot(Result[:,0],Result[:,3],linewidth=3, label='torsional stiffness',color='blue',marker='o',markersize=8,linestyle='dashed')
bend, = ax4.plot(Result[:,0],Result[:,4],linewidth=3, label='bending stiffness',color='red',marker='^',markersize=8,linestyle='dashed')
ax4.set_ylabel('stiffness (N.m²)',fontsize=20)
ax4.tick_params(labelsize=20)
plt.legend(handles=[div,tors,bend],loc=6, prop={'size': 25})

plt.grid()
plt.show()

