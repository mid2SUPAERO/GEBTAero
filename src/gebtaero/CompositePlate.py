# coding=UTF-8
import numpy as np
import math
from .utils import *
from .OrthoMaterial import *
from .IsoMaterial import *

## Class interfacing the solver with 3D FEM calculix computation to obtain the cross section parameter from a composite plate
class CompositePlate:
    """
    Class interfacing the solver with 3D FEM calculix computation to obtain stiffness matrix
    from a composite plate
    """        
    def __init__(self,Chord=1,OffsetY=0.,OffsetZ=0.):
        ## Width of the plate
        self.Chord = Chord  
        ## A list containing some CompositePly
        self.Layup = []
        ## A list containing the differents materials used in the plata (either IsoMaterial or OrthoMaterial)
        self.Materials = []
        ## A list containing the differents fiber orientations of the CompositePly of the plate
        self.Orientations = []
        ## The total thickness of the plate
        self.TotThickness = 0.
        ## Y coordinate of the plate center in the Frame B
        self.OffsetY = OffsetY
        ## Y coordinate of the plate center in the Frame B
        self.OffsetZ = OffsetZ
        
    ## Add a CompositePly to the plate    
    def AppendPly(self,Ply):
        self.Layup.append(Ply)
        if(Ply.GetMaterial() not in self.Materials):
            self.Materials.append(Ply.GetMaterial())
        if(Ply.GetOrientation() not in self.Orientations):
            self.Orientations.append(Ply.GetOrientation())
        self.TotThickness = self.TotThickness + Ply.GetThickness()
        
    def GetLayup(self):
        return self.Layup
        
    def GetTotThickness(self):
        return self.TotThickness 
    
    def GetOffsets(self):
        return [self.OffsetY,self.OffsetZ]    
    
    ## Create the input file for cgx preprocessor  
    #@param TypeElem the type of finite element used for the computation. Supported : he20r, he20, he8i, he8, pe15 (see ccx doc)
    #@param NbElemX the number of finite element across the beam direction (for a constant cross section, 1 element is enough)
    #@param NBElemY the number of finite element across the width
    #@param NBElemPly the number of finite element in a composite ply
    def CreateFbdFile(self,TypeElem,NbElemX,NbElemY,NbElemPly):
        NbPly = len(self.GetLayup())
        if (NbPly <1):
            raise RuntimeError("A composite plate must have at lest one composite ply")
        if(TypeElem in ["he20r","he20"]):
            divX = 2*NbElemX
            divY = 2*NbElemY
            divPly = 2*NbElemPly
        elif(TypeElem in ["he8i","he8","pe15"]):   
            divX = NbElemX
            divY = NbElemY
            divPly = NbElemPly
        else:
            raise RuntimeError("The type of element for the plate homogeneisation is not recognize")
        # determination of RVE dimensions
        NbPly = len(self.GetLayup())
        # ~ Lx = round(divX*0.5*(self.Chord/divY+self.TotThickness/(NbPly*divPly)),10)
        # ~ Ly = round(self.Chord,10)
        # ~ Lz = round(self.GetTotThickness(),10)
        Lx = round(divX*0.5*(1.0/divY+1.0*self.TotThickness/self.Chord/(NbPly*divPly)),16)
        Ly = round(1.0,16)
        Lz = round(self.GetTotThickness()/self.Chord,16)
        file=open("Plate.fbd","w")
        file.write("PNT ! {0} {1} {2}\n".format(-0.5*Lx,-0.5*Ly+self.OffsetY/self.Chord,-0.5*Lz+self.OffsetZ/self.Chord))
        file.write("SWEP all all tra {0} 0. 0. {1}\n".format(Lx,divX))
        file.write("SWEP all all tra 0. {0} 0. {1}\n".format(Ly,divY))
        # ~ file.write("SWEP all plan{0} tra 0. 0. {1} {2}\n".format(1,round(self.GetLayup()[0].GetThickness(),10),divPly))
        file.write("SWEP all plan{0} tra 0. 0. {1} {2}\n".format(1,round(self.GetLayup()[0].GetThickness()/self.Chord,16),divPly))
        file.write("SETA ply1 b all\n")
        for i in range (NbPly-1):
            file.write("SETA plym b all\n")            
            # ~ file.write("SWEP plan{0} plan{1} tra 0. 0. {2} {3}\n".format(i+1,i+2,round(self.GetLayup()[i+1].GetThickness(),10),divPly))
            file.write("SWEP plan{0} plan{1} tra 0. 0. {2} {3}\n".format(i+1,i+2,round(self.GetLayup()[i+1].GetThickness()/self.Chord,16),divPly))
            # ~ file.write("SWEP plan{0} plan{1} tra 0. 0. {2} {3}\n".format(i+1,i+2,round(self.GetLayup()[i+1].GetThickness(),10),divPly))
            # ~ file.write("SWEP plan{0} plan{1} tra 0. 0. {2} {3}\n".format(i+1,i+2,round(self.GetLayup()[i+1].GetThickness()/self.Chord,16),divPly))
            file.write("SETA ply{0} b all\n".format(i+2))
            file.write("SETR ply{0} b plym\n".format(i+2))
        file.write('MERG p all\n')    
        file.write('MERG l all\n')    
        file.write('MERG s all\n') 
        file.write("ELTY all "+TypeElem+"\n")
        file.write("MESH all\n")
        file.write("SEND all abq\n")
        for i in range (NbPly):
            file.write("SEND ply{0} abq nam\n".format(i+1))   
        file.close
        
    ## Create the input file for ccx solver
    #@param Stress True = output the stress tensor for each elementary load case, False = no stress output
    #@param PlaneSection True = warping of the cross section is not allowed, False = warping of the cross section is allowed
    #@param Disp  Determine the set of elementary load cases : 0=all the 4 load cases, 1 = traction, 2 = warping, 3 = bending span-wise, 4 = bending chord-wise.        
    def CreateInpFile(self,Stress=False,PlaneSection=False,Disp=0):
        NbPly = len(self.Layup)
        NbMat = len(self.Materials)
        NbOri = len(self.Orientations)
        if (NbPly <1):
            raise RuntimeError("A composite plate must have at lest one composite ply")
        file=open("Plate.inp","w")
        file.write("*include,input=all.msh\n")
        #~ file.write("*include,input=all.nam\n")
        for i in range(NbPly):
            file.write("*include,input=ply"+str(i+1)+".nam\n")
        file.write("*include,input=periodic.equ\n")    
        # materials
        for i in range(NbMat):
            Mat = self.Materials[i]
            if isinstance(Mat,OrthoMaterial):
                El = Mat.GetOrtho()[0]
                Et = Mat.GetOrtho()[1]
                Nult = Mat.GetOrtho()[2]
                Nutl = round(float(Et)/El*Nult,6)
                Glt = Mat.GetOrtho()[3]
                Gtl = round(float(Et)/(2+2*Nutl),0)
                file.write("*material,name=MAT"+str(i+1)+"\n")
                file.write("*elastic,type=ENGINEERING CONSTANTS\n")
                file.write(str(El)+","+str(Et)+","+str(Et)+","+str(Nult)+","+str(Nutl)+","+str(Nutl)+","+str(Glt)+","+str(Glt)+"\n"+str(Gtl)+"\n")
            elif isinstance(Mat,IsoMaterial):
                E = Mat.GetIso()[0]
                Nu = Mat.GetIso()[1]
                file.write("*material,name=MAT"+str(i+1)+"\n")
                file.write("*elastic,type=ISO\n")
                file.write(str(E)+","+str(Nu)+'\n')
        # orientation
        for i in range(NbOri):
            alpha = np.pi/180*self.Orientations[i]
            cos = str(round(math.cos(alpha),6))
            sin = str(round(math.sin(alpha),6))
            msin = str(round(-math.sin(alpha),6))
            file.write("*orientation,name=Ori"+str(i+1)+"\n")
            file.write(cos+","+sin+",0.,"+msin+","+cos+",0.\n")
        # elements affectation
        for i in range(NbPly):
            iMat = self.Materials.index(self.Layup[i].GetMaterial())+1
            iOri = self.Orientations.index(self.Layup[i].GetOrientation())+1
            if isinstance(self.Materials[iMat-1],OrthoMaterial):
                file.write("*solid section, material=MAT"+str(iMat)+",elset=Eply"+str(i+1)+",orientation=Ori"+str(iOri)+"\n")
            elif isinstance(self.Materials[iMat-1],IsoMaterial):
                file.write("*solid section, material=MAT"+str(iMat)+",elset=Eply"+str(i+1)+"\n")     
        file.write("\n")
        if PlaneSection:
            file.write("*boundary\nnxl,1,1\n\n")
        #~ # elementary beam stress loading for flexibility matrix computation
        if Disp == 0:
            if Stress:    
                file.write("*step\n*static\n*cload,OP=NEW\nnstrain,1,1.\n*node print,NSET=nstrain\nU\n*node print,NSET=ncurv\nU\n*node file,NSET=nxl\nS\n*end step\n\n")
                file.write("*step\n*static\n*cload,OP=NEW\nncurv,1,1.\n*node print,NSET=nstrain\nU\n*node print,NSET=ncurv\nU\n*node file,NSET=nxl\nS\n*end step\n\n")
                file.write("*step\n*static\n*cload,OP=NEW\nncurv,2,1.\n*node print,NSET=nstrain\nU\n*node print,NSET=ncurv\nU\n*node file,NSET=nxl\nS\n*end step\n\n")
                file.write("*step\n*static\n*cload,OP=NEW\nncurv,3,1.\n*node print,NSET=nstrain\nU\n*node print,NSET=ncurv\nU\n*node file,NSET=nxl\nS\n*end step\n\n")
            else:
                file.write("*step\n*static\n*cload,OP=NEW\nnstrain,1,1.\n*node print,NSET=nstrain\nU\n*node print,NSET=ncurv\nU\n*end step\n\n")
                file.write("*step\n*static\n*cload,OP=NEW\nncurv,1,1.\n*node print,NSET=nstrain\nU\n*node print,NSET=ncurv\nU\n*end step\n\n")
                file.write("*step\n*static\n*cload,OP=NEW\nncurv,2,1.\n*node print,NSET=nstrain\nU\n*node print,NSET=ncurv\nU\n*end step\n\n")
                file.write("*step\n*static\n*cload,OP=NEW\nncurv,3,1.\n*node print,NSET=nstrain\nU\n*node print,NSET=ncurv\nU\n*end step\n\n")
        elif Disp == 1:
            file.write("*step\n*static\n*cload,OP=NEW\nnstrain,1,1.\n*node file,NSET=Nall\nU\n*end step\n\n")
        elif Disp == 2:
            file.write("*step\n*static\n*cload,OP=NEW\nncurv,1,1.\n*node file,NSET=Nall\nU\n*end step\n\n")
        elif Disp == 3:
            file.write("*step\n*static\n*cload,OP=NEW\nncurv,2,2.\n*node file,NSET=Nall\nU\n*end step\n\n")
        elif Disp == 4:
            file.write("*step\n*static\n*cload,OP=NEW\nncurv,3,3.\n*node file,NSET=Nall\nU\n*end step\n\n")
        file.close()
        
    def CreateVtkFile(self,ChordScale=False):
        L = self.Chord
        h = self.TotThickness
        p0x = -0.5*L+self.OffsetY
        p0y = -0.5*h+self.OffsetZ
        p1x = -0.5*L+self.OffsetY
        p1y = 0.5*h+self.OffsetZ
        p2x = 0.5*L+self.OffsetY
        p2y = 0.5*h+self.OffsetZ
        p3x = 0.5*L+self.OffsetY
        p3y = -0.5*h+self.OffsetZ
        
        if ChordScale:
            p0x = p0x/L
            p0y = p0y/L
            p1x = p1x/L
            p1y = p1y/L
            p2x = p2x/L
            p2y = p2y/L
            p3x = p3x/L
            p3y = p3y/L

        file=open("Plate.vtk","w")
        file.write("# vtk DataFile Version 4.2\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS 4 float\n")
        file.write("{0} {1} 0\n".format(p0x,p0y))
        file.write("{0} {1} 0\n".format(p1x,p1y))
        file.write("{0} {1} 0\n".format(p2x,p2y))
        file.write("{0} {1} 0\n".format(p3x,p3y))
        file.write("\nPOLYGONS 1 5\n4 0 1 2 3\n")
     
    ## Compute the CrossSection::MassMatrix analytically 
    #@return Mu Mass per unit length (kg/m)
    #@return I22 mass moment of inertia around yB (bending span-wise, kg.m) 
    #@return I33 mass moment of inertia around zB (bending chord-wise, kg.m) 
    #@return I23 product of inertia in plane (yB,zB)
    #@return Ycg Y coordinate of the Center of Gravity in frame B
    #@return Zcg Z coordinate of the Center of Gravity in frame B
    def ComputeMassMatrix(self,OffsetY=0.,OffsetZ=0.):
        Mu = 0.
        I22 = 0.
        I33 = 0.
        Zcg = 0.
        Iply = []
        ZLayup = [0.]
        NbPly = len(self.Layup)
        Chord = self.Chord
        Lz = self.TotThickness
        for i in range(NbPly):
            Rho= self.Layup[i].GetMaterial().GetDensity()
            if Rho==None:
                raise RuntimeError('You must define a density before calculting the Mass Matrix')
            Ep = self.Layup[i].GetThickness()
            MuPly = Chord*Ep*Rho
            Mu = Mu + MuPly
            ZLayup.append(ZLayup[-1]+Ep)
            ZplyCenter = round(0.5*(ZLayup[-1]+ZLayup[-2])-0.5*Lz-OffsetZ,16)
            I22 = I22 + MuPly/12*Ep**2+MuPly*ZplyCenter**2
            Zcg = Zcg + MuPly*ZplyCenter
        I33 = Mu/12*Chord**2 + Mu*OffsetY**2
        I23 = Mu*OffsetY*OffsetZ
        Zcg = round(1+Zcg/Mu,12)-1
        Ycg = -OffsetY
        return [Mu,I22,I33,I23,Ycg,Zcg]
        
    ##  Interface with the function #utils::CreatePeriodicEq
    def CreatePeriodicEq(self):
        [self.nstrain,self.ncurv,self.x0,self.Lx,self.nodes,self.elements]=CreatePeriodicEq("all.msh")

    ## Launch the cgx postprocessor with a particular elementary load case.
    #@param TypeElem the type of finite element used for the computation. Supported : he20r, he20, he8i, he8, pe15
    #@param NbElemX the number of finite element in the beam direction (for a constant cross section, 1 element is enough)
    #@param NBElemY the number of finite element across the width
    #@param NBElemPly the number of finite element in a composite ply
    #@param DefType  Determine the set of elementary load cases : 0=all the 4 load cases, 1 = traction, 2 = warping, 3 = bending span-wise, 4 = bending chord-wise.
    #@param PlaneSection : True = warping of the cross section is not allowed, False = warping of the cross section is allowed        
    def DisplaySectionDeformation(self,TypeElem,NbElemX,NbElemY,NbElemPly,DefType,PlaneSection=False):
        self.CreateFbdFile(TypeElem,NbElemX,NbElemY,NbElemPly)
        RunFbdFile("Plate.fbd")
        self.CreatePeriodicEq()
        self.CreateInpFile(Disp=DefType,PlaneSection=PlaneSection)
        Status = RunInpFile("Plate")[0]
        if Status !=0:
            raise RuntimeError("calculix computation failure")
        RunFrdFile("Plate.frd")
