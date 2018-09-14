# coding=UTF-8
import numpy as np
from .utils import *
from .IsoMaterial import *

## This class is used to determine the Mass matrix and the flexibility matrix of a cross section. It could be done analytically or with the 3D FEM solver ccx.
class CrossSection:
    
    ## The constructor initialize both matrices to zero.
    def __init__(self):
        ## Mass matrix of the cross section (see GEBTAero paper)
        self.MassMatrix = np.zeros([6,6])
        ## Flexibility matrix of the cross section (see GEBTAero paper)
        self.FlexibilityMatrix = np.zeros([6,6])
        ## Coordinate of the elastic center in the plane (yB,Zb)
        self.ElasticCenter = None
        
    ## Set the Mass matrix using classical isotropic test cases parameters
    #@param Mu Mass per unit length (kg/m)
    #@param I22 mass moment of inertia around yB (bending span-wise, kg.m) 
    #@param I33 mass moment of inertia around zB (bending chord-wise, kg.m) 
    #@param I23 product of inertia in plane (yB,zB)
    #@param Nu : Distance between moment calculation point and center of gravity. Positif is the CG is near leading edge
    #@param Zcg : coordinate Z of the center of gravity in frame B
    def SetMassMatrix(self, Mu, I22, I33, I23, Nu, Zcg=0.):
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

    ## Interface to CompositePlate::ComputeMassMatrix
    def SetMassMatrixByPlate(self,Plate):
        [OffsetY,OffsetZ] = Plate.GetOffsets()
        [Mu,I22,I33,I23,Ycg,Zcg]=Plate.ComputeMassMatrix(OffsetY,OffsetZ)
        self.SetMassMatrix(Mu,I22,I33,I23,-Ycg,Zcg=Zcg)   
    
    ## Interface to CompositeBox::ComputeMassMatrix
    def SetMassMatrixByBox(self,Box):
        [OffsetY,OffsetZ] = Box.GetOffsets()
        [Mu,I22,I33,I23,Ycg,Zcg]=Box.ComputeMassMatrix(OffsetY,OffsetZ)
        self.SetMassMatrix(Mu,I22,I33,I23,-Ycg,Zcg=Zcg)   
    
    ## Analytical definition of the CrossSection::MassMatrix in case of a rectangular section isotropic beam
    #@param h Section height
    #@param L Section width
    #@param Mat Material of the beam (must be IsoMaterial)
    #@param NeglectI22 True = the value of I22 is forced to zero, False = nothing appends
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
            
    ## Interface to Mesh::ComputeMassMatrixFromMesh
    #@param Mesh The Mesh used for the calculation (inp format)
    #@param AtElasticCenter Compute the Elastic center and then compute the mass matrix at the elastic center
    #@param SymY True = force the horizontal symetry of the cross section
    #@param ChordScale True = the Mesh have a dimension 1 in the chord direction (yB) the mass matrix obtained is scaled afterward with the parameter Mesh::Chord
    #@param verbosity 0 = no log, 1 = output the log of the computation
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
        
    ## Analytical definition of the CrossSection::FlexibilityMatrix for a rectangular section isotropic beam    
    #@param h Section height
    #@param L Section width
    #@param Mat Material of the beam (must be IsoMaterial)
    #@param RigidEIg3 True = supress the bending ddl in the lag axis
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

    ## Set the CrossSection::FlexibilityMatrix using isotropic beam stiffness values.
    #@param EIg2 bending stiffness span-wise
    #@param EIg3 bending stiffness chord-wise
    #@param GJ torsional stiffness
    def SetFlexibilityMatrixByIsotropicValues(self, EIg2, EIg3, GJ):
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

    ## Set the CrossSection::FlexibilityMatrix using a periodic 3DFEM calculation with a CompositePlate
    #@param Plate the CompositePlate used for the calculation
    #@param TypeElem the type of finite element used for the computation. Supported : he20r, he20, he8i, he8, pe15 (see ccx doc)
    #@param NbElemX the number of finite element across the beam direction (for a constant cross section, 1 element is enough)
    #@param NBElemY the number of finite element across the width
    #@param NBElemPly the number of finite element in a composite ply    
    #@param RigidX True = supress the traction ddl of the beam
    #@param RigidZ True = supress the bending in the lag axis ddl
    def SetFlexibilityMatrixByPlate(self,Plate,TypeElem,NbElemX,NbElemY,NbElemPly,RigidX=False,RigidZ=False):
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
        
    ## Set the CrossSection::FlexibilityMatrix using a periodic 3DFEM calculation with a CompositeBox
    #@param Box the CompositeBox used for the calculation
    #@param TypeElem the type of finite element used for the computation. Supported : he20r, he20, he8i, he8, pe15 (see ccx doc)
    #@param NbElemX the number of finite element across the beam direction (for a constant cross section, 1 element is enough)
    #@param NBElemYZ the number of finite element along the wall directions
    #@param NBElemPly the number of finite element in a composite ply    
    #@param RigidX True = supress the traction ddl of the beam
    #@param RigidZ True = supress the bending in the lag axis ddl        
    def SetFlexibilityMatrixByBox(self,Box,TypeElem,NbElemX,NbElemYZ,NbElemPly,RigidX=False,RigidZ=False):
        Box.CreateFbdFile(TypeElem,NbElemX,NbElemYZ,NbElemPly)
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
        
    ## Set the CrossSection::FlexibilityMatrix using Mesh
    #@param Mesh The Mesh used for the calculation (inp format)
    #@param PlaneSection True = warping of the cross section is not allowed, False = warping of the cross section is allowed
    #@param AtElasticCenter Compute the Elastic center and then compute the mass matrix at the elastic center
    #@param RigidX True = supress the traction ddl of the beam
    #@param RigidZ True = supress the bending in the lag axis ddl
    #@param ChordScale True = the Mesh have a dimension 1 in the chord direction (yB) the mass matrix obtained is scaled afterward with the parameter Mesh::Chord      
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
