# coding=UTF-8
import numpy as np
from .utils import *
from .IsoMaterial import *

class CrossSection:
    """
    class containing the flexibility and mass matrix of a wing section
    Mass matrix is defined analytically, flexibility matrix is defined either analytically 
    or using 3D FEM solver Calculix
    """
    def __init__(self):
        self.MassMatrix = np.zeros([6,6])
        self.FlexibilityMatrix = np.zeros([6,6])
        self.ElasticCenter = None
        
    def SetMassMatrix(self, Mu, I22, I33, I23, Nu, Zcg=0.):
        """
        @param float Mu : Mass per unit length
        @param float I11 : Mass moment of inertia about wingspan axis
        @param float Nu : Distance between moment calculation point and center of gravity. Positif is the CG is near leading edge
        """
        if (Mu < 0):
            raise ValueError ('Mu must be positive')
        else :
            self.MassMatrix[0][0] = self.MassMatrix[1][1] = self.MassMatrix[2][2] = Mu
        if (I22 < 0 or I33 < 0):
            raise ValueError ('I11 must be positive')
        else :
            Ycg = -Nu
            self.MassMatrix[3][3] = I22+I33
            #Inertia around airfoil axis is negligible but cannot be zero (inversion of mass matrix in the solver)
            self.MassMatrix[4][4] = I22
            self.MassMatrix[5][5] = I33
            self.MassMatrix[4][5] = self.MassMatrix[5][4] = -I23
            self.MassMatrix[0][5] = self.MassMatrix[5][0] = -Mu*Ycg
            self.MassMatrix[2][3] = self.MassMatrix[3][2]  = Mu*Ycg
            self.MassMatrix[0][4] = self.MassMatrix[4][0] = Mu*Zcg
            self.MassMatrix[1][3] = self.MassMatrix[3][1]  = -Mu*Zcg

    def SetMassMatrixByPlate(self,Plate):
        [OffsetY,OffsetZ] = Plate.GetOffsets()
        [Mu,I22,I33,I23,Ycg,Zcg]=Plate.ComputeMassMatrix(OffsetY,OffsetZ)
        self.SetMassMatrix(Mu,I22,I33,I23,-Ycg,Zcg=Zcg)   
    
    def SetMassMatrixByBox(self,Box):
        [OffsetY,OffsetZ] = Box.GetOffsets()
        [Mu,I22,I33,I23,Ycg,Zcg]=Box.ComputeMassMatrix(OffsetY,OffsetZ)
        self.SetMassMatrix(Mu,I22,I33,I23,-Ycg,Zcg=Zcg)   
    
    def SetMassMatrixByRectBeamValues(self,h,L,Mat,NeglectI22=False):
        if not isinstance(Mat,IsoMaterial):
            raise TypeError("the material must be isotropic") 
        else:       
            Rho = Mat.GetDensity()
            I33 = Rho*L**3*h/12.   
            if NeglectI22:
                I22=1e-4*I33
            else:
                I22 = Rho*L*h**3/12.
             
            Mu=Rho*L*h
            self.SetMassMatrix(Mu,I22,I33,0.,0.)
            
        
    def SetMassMatrixByMesh(self,Mesh,AtElasticCenter=False,SymY=False,ChordScale=False,verbosity=0):
        if AtElasticCenter:
            if self.ElasticCenter is None:
                [nstrain,ncurv,x0,Lx,nodes,elements]=CreatePeriodicEq(Mesh.GetMeshFile())
                Mesh.CreateInpFile(Stress=False,PlaneSection=False)
                Status = RunInpFile("Ext")[0]
                if Status !=0:
                    raise RuntimeError("calculix computation failure")
                self.FlexibilityMatrix = ReadFlexibilityFromDisp("Ext.dat",nstrain,ncurv,Lx,1e-6)
                RemoveFiles()
            self.ElasticCenter = ComputeElasticCenterFromFlexMat(self.FlexibilityMatrix) 
            #backup Offsets
            OY = Mesh.OffsetY    
            OZ = Mesh.OffsetZ
            Mesh.OffsetY = self.ElasticCenter[0]
            Mesh.OffsetZ = self.ElasticCenter[1]
        [Mu,I22,I33,I23,Ycg,Zcg]=Mesh.ComputeMassMatrixFromMesh(Mesh.nodes,Mesh.x0,Mesh.elements)
        if verbosity ==1:
            print("Mu = {0}, I22 = {1}, I33 = {2}, I23 = {3}, Ycg = {4}, Zcg = {5}".format(Mu,I22,I33,I23,Ycg,Zcg))
        #common part
        if SymY:
            self.SetMassMatrix(Mu,I22,I33,I23,-Ycg) 
        else:    
            self.SetMassMatrix(Mu,I22,I33,I23,-Ycg,Zcg)     
        if AtElasticCenter:    
            #Offsets restauration
            Mesh.OffsetY = OY
            Mesh.OffsetZ = OZ
        if ChordScale:
            # case of adimensional chord calculation
            self.MassMatrix[0:3,0:3] = self.MassMatrix[0:3,0:3] *(Mesh.Chord**2)
            self.MassMatrix[3:6,3:6] = self.MassMatrix[3:6,3:6] *(Mesh.Chord**4)
            self.MassMatrix[0:3,3:6] = self.MassMatrix[0:3,3:6] *(Mesh.Chord**3)
            self.MassMatrix[3:6,0:3] = self.MassMatrix[3:6,0:3] *(Mesh.Chord**3)
                
            
            
            
    
    def GetMassMatrix(self):
        return self.MassMatrix
        
    def SetFlexibilityMatrixByRectBeamValues(self,h,L,Mat,RigidEIg3=False):
        if isinstance(Mat,IsoMaterial):
            E = Mat.GetIso()[0]
            Nu = Mat.GetIso()[1]
            G = E/(2*(1+Nu))
            Ig2 = L*h**3/12.
            if RigidEIg3:
                Ig3 = 0
            else:
                Ig3 = L**3*h/12.
            # Roark formula p 406    
            a = 0.5*L
            b = 0.5*h
            J = a*b**3*(16./3-3.36*b/a*(1-b**4/(12*a**4)))
        else:
            raise TypeError("the material must be isotropic")    
        self.SetFlexibilityMatrixByIsotropicValues(E*Ig2,E*Ig3,G*J)

    def SetFlexibilityMatrixByIsotropicValues(self, EIg2, EIg3, GJ):
        """
        Set the flexibility matrix using isotropic vlaues namely the bending stiffness
        spanwise EIg2, the bending stiffness chordwise EIg3 and the torsionnal stiffness
        GJ.
        @param float EIg2 : bending stiffness spanwise
        @param float EIg3 : bending stiffness chordwise
        @param float GJ : torsionnal stiffness
        """
        if (EIg2 <0):
            raise ValueError ('bending stiffness must be positive')
        elif (EIg2 > 0):
         self.FlexibilityMatrix[4][4] = 1./EIg2

        if (EIg3 <0):
            raise ValueError ('bending stiffness must be positive')
        elif (EIg3 > 0) :
            self.FlexibilityMatrix[5][5]= 1./EIg3		

        if (GJ <0) :
            raise ValueError ('bending stiffness must be positive')
        elif (GJ > 0):
            self.FlexibilityMatrix[3][3] = 1./GJ				

    def SetFlexibilityMatrixByPlate(self,Plate,TypeElem,NbElemX,NbElemY,NbElemPly,RigidX=False,RigidZ=False):
        """
        Set the flexibility matrix using homogeneisation process.
        First a 4*4 Stiffness matrix is computed, then inverted and coefficients are
        placed in the final 6*6 Flexibility Matrix
        """
        Plate.CreateFbdFile(TypeElem,NbElemX,NbElemY,NbElemPly)
        RunFbdFile("Plate.fbd")
        Plate.CreatePeriodicEq()
        Plate.CreateInpFile()
        Status = RunInpFile("Plate")[0]
        if Status !=0:
            raise RuntimeError("calculix computation failure")            
        self.FlexibilityMatrix=ReadFlexibilityFromDisp("Plate.dat",Plate.nstrain,Plate.ncurv,Plate.Lx,1e-10 ,RigidX,RigidZ)
        # case of adimensional chord calculation
        self.FlexibilityMatrix[0:3,0:3] = self.FlexibilityMatrix[0:3,0:3] *(Plate.Chord**-2)
        self.FlexibilityMatrix[3:6,3:6] = self.FlexibilityMatrix[3:6,3:6] *(Plate.Chord**-4)
        self.FlexibilityMatrix[0:3,3:6] = self.FlexibilityMatrix[0:3,3:6] *(Plate.Chord**-3)
        self.FlexibilityMatrix[3:6,0:3] = self.FlexibilityMatrix[3:6,0:3] *(Plate.Chord**-3)
        RemoveFiles()
        
    def SetFlexibilityMatrixByBox(self,Box,TypeElem,NbElemX,NbElemY,NbElemPly,RigidX=False,RigidZ=False):
        """
        Set the flexibility matrix using homogeneisation process.
        First a 4*4 Stiffness matrix is computed, then coefficients are
        placed in the final 6*6 Flexibility Matrix and scaled using the box width
        """
        Box.CreateFbdFile(TypeElem,NbElemX,NbElemY,NbElemPly)
        RunFbdFile("Box.fbd")
        Box.CreatePeriodicEq()
        Box.CreateInpFile()
        Status = RunInpFile("Box")[0]
        if Status !=0:
            raise RuntimeError("calculix computation failure")            
        self.FlexibilityMatrix=ReadFlexibilityFromDisp("Box.dat",Box.nstrain,Box.ncurv,Box.Lx,1e-10 ,RigidX,RigidZ)
        # case of adimensional chord calculation
        self.FlexibilityMatrix[0:3,0:3] = self.FlexibilityMatrix[0:3,0:3] *(Box.Width**-2)
        self.FlexibilityMatrix[3:6,3:6] = self.FlexibilityMatrix[3:6,3:6] *(Box.Width**-4)
        self.FlexibilityMatrix[0:3,3:6] = self.FlexibilityMatrix[0:3,3:6] *(Box.Width**-3)
        self.FlexibilityMatrix[3:6,0:3] = self.FlexibilityMatrix[3:6,0:3] *(Box.Width**-3)
        RemoveFiles()
        
    def SetFlexibilityMatrixByMesh(self,Mesh,PlaneSection=False,AtElasticCenter=False,RigidX=False,RigidZ=False,ChordScale=False):
        Mesh.CreatePeriodicEq()
        Mesh.CreateInpFile(Stress=False,PlaneSection=PlaneSection)
        Status = RunInpFile("Ext")[0]
        if Status !=0:
            raise RuntimeError("calculix computation failure")            
        self.FlexibilityMatrix = ReadFlexibilityFromDisp("Ext.dat",Mesh.nstrain,Mesh.ncurv,Mesh.Lx,1e-6,RigidX,RigidZ)
        RemoveFiles()
        if AtElasticCenter:
            if self.ElasticCenter is None:
                self.ElasticCenter = ComputeElasticCenterFromFlexMat(self.FlexibilityMatrix) 
            [nstrain,ncurv,x0,Lx,nodes,elements]=CreatePeriodicEq(Mesh.GetMeshFile(),OffsetY=self.ElasticCenter[0],OffsetZ=self.ElasticCenter[1])
            Mesh.CreateInpFile(Stress=False,PlaneSection=PlaneSection)
            Status = RunInpFile("Ext")[0]
            if Status !=0:
                raise RuntimeError("calculix computation failure")                        
            self.FlexibilityMatrix = ReadFlexibilityFromDisp("Ext.dat",nstrain,ncurv,Lx,1e-6,RigidX,RigidZ)
            RemoveFiles()
        if ChordScale:
            # case of adimensional chord calculation
            self.FlexibilityMatrix[0:3,0:3] = self.FlexibilityMatrix[0:3,0:3] *(Mesh.Chord**-2)
            self.FlexibilityMatrix[3:6,3:6] = self.FlexibilityMatrix[3:6,3:6] *(Mesh.Chord**-4)
            self.FlexibilityMatrix[0:3,3:6] = self.FlexibilityMatrix[0:3,3:6] *(Mesh.Chord**-3)
            self.FlexibilityMatrix[3:6,0:3] = self.FlexibilityMatrix[3:6,0:3] *(Mesh.Chord**-3)
            
        
    def GetFlexibilityMatrix(self):
        return self.FlexibilityMatrix
