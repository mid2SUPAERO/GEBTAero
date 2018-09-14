# coding=UTF-8
import numpy as np
import math
from .utils import *
from .OrthoMaterial import *
from .IsoMaterial import *

## Class interfacing the solver with 3D FEM calculix computation to obtain the cross section parameter from a composite box
class CompositeBox:
    ## Composite box constructor; It contains 4 CompositePlate; modify the chord of left and right wall  
    def __init__(self,Left,Right,Up,Down,Width,Height,OffsetY=0.,OffsetZ=0.):
        ## Left wall of the box (CompositePlate)
        self.Left = Left    
        ## Right wall of the box (CompositePlate)
        self.Right = Right
        ## Up wall of the box (CompositePlate)
        self.Up = Up
        ## Down wall of the box (CompositePlate)
        self.Down = Down
        ## Width of the box (m)
        self.Width = Width
        ## Height height of the box (m)
        self.Height = Height
        ## Y coordinate of the box center in the Frame B
        self.OffsetY = OffsetY
        ## Z coordinate of the box center in the Frame B
        self.OffsetZ = OffsetZ
        self.Up.Chord = Width
        self.Down.Chord = Width
        self.Left.Chord = Height - self.Up.GetTotThickness()- self.Down.GetTotThickness()
        self.Right.Chord = Height - self.Up.GetTotThickness()- self.Down.GetTotThickness()
        
    def GetOffsets(self):
        return [self.OffsetY,self.OffsetZ]    
        
    def GetWidth(self):
        return self.Width
        
    def GetHeight(self):
        return self.Height
        
    ## Create the input file for cgx preprocessor  
    #@param TypeElem the type of finite element used for the computation. Supported : he20r, he20, he8i, he8, pe15 (see ccx doc)
    #@param NbElemX the number of finite element in the beam direction (for a constant cross section, 1 element is enough)
    #@param NBElemYZ the number of finite element along the wall directions
    #@param NBElemPly the number of finite element in a composite ply  
    def CreateFbdFile(self,TypeElem,NbElemX,NbElemYZ,NbElemPly):
        if(TypeElem in ["he20r","he20"]):
            divX = 2*NbElemX
            divYZ = 2*NbElemYZ
            divPly = 2*NbElemPly
        elif(TypeElem in ["he8i","he8","pe15"]):   
            divX = NbElemX
            divYZ = NbElemYZ
            divPly = NbElemPly
        else:
            raise RuntimeError("The type of element for the plate homogeneisation is not recognize")
        # determination of RVE dimensions
        if min(len(self.Left.Layup),len(self.Right.Layup),len(self.Up.Layup),len(self.Down.Layup)) ==0:
            raise RuntimeError("One of the box sides has no ply defined")
        TL = self.Left.Layup[0].Thickness
        TR = self.Right.Layup[0].Thickness
        TU = self.Up.Layup[0].Thickness
        TD = self.Down.Layup[0].Thickness
        H = self.Height
        W = self.Width
        n=[]
        n.append(len(self.Up.Layup))
        n.append(len(self.Left.Layup))
        n.append(len(self.Down.Layup))
        n.append(len(self.Right.Layup))
        if min(TL,TR,TU,TD)>0:
            pass
        else:
            raise RuntimeError("One of the box sides has a negative or zero thickness")    
        Lx = round(divX*0.25*(TL+TR+TU+TD)/NbElemPly/self.Width,16)
        Ly = 1.0
        Lz = round(H/W,16)
        file=open("Box.fbd","w")
        # creation of the left up reference point
        file.write("SETO up\nSETO up0\nSETO left\n")
        file.write("PNT p1 {0} {1} {2}\n".format(-0.5*Lx,0.5+self.OffsetY/W,0.5*Lz-self.Up.GetTotThickness()/W+self.OffsetZ/W))
        file.write("SETC left\n")
        #creation of the upper side
        # dividing the upper side in accordance to the ply division of left and right side
        
        
        for i in range(n[1]): #left side
            file.write("SWEP up{2} up{3} tra 0 {0} 0 {1}\n".format(-round(self.Left.Layup[n[1]-1-i].GetThickness()/W,16),divPly,i,i+1))
        file.write("SWEP up{2} up{3} tra 0 {0} 0 {1}\n".format(round(-1+(self.Left.GetTotThickness()+self.Right.GetTotThickness())/W,16),divYZ,n[1],n[1]+1))
        for i in range(n[3]): #right side
            file.write("SWEP up{2} up{3} tra 0 {0} 0 {1}\n".format(-round(self.Right.Layup[i].GetThickness()/W,16),divPly,n[1]+1+i,n[1]+2+i))
        
        
        file.write("SETO ply1up\n")
        file.write("SWEP up plan1up tra 0 0 {0} {1}\n".format(self.Up.Layup[0].GetThickness()/W,divPly))
        file.write("SWEP ply1up ply1up tra {0} 0 0 {1}\n".format(Lx,divX))
        file.write("SETR ply1up p all\nSETR ply1up l all\nSETR ply1up s all\n")
        file.write("SETC ply1up\n")        
        for i in range (n[0]-1):
            file.write("SETO ply{0}up\n".format(i+2))
            file.write("SWEP plan{0}up plan{1}up tra 0 0 {2} {3}\n".format(i+1,i+2,round(self.Up.Layup[i+1].GetThickness()/W,16),divPly))
            file.write("SWEP ply{0}up ply{1}up tra {2} 0 0 {3}\n".format(i+2,i+2,Lx,divX))
            file.write("SETR ply{0}up p all\nSETR ply{0}up l all\nSETR ply{0}up s all\n".format(i+2))
            file.write("SETC ply{0}up\n".format(i+2))        
        #creation of the left side
        file.write("SETC up\nSETO left\n") 
        file.write("SWEP left left tra 0 0 {0} {1}\n".format(round((-H+self.Up.GetTotThickness()+self.Down.GetTotThickness())/W,16),divYZ))
        file.write("SETO ply{0}left\n".format(n[1]))
        file.write("SWEP left plan{2}left tra 0 {0} 0 {1}\n".format(round(-self.Left.Layup[n[1]-1].GetThickness()/W,16),divPly,n[1]))
        file.write("SWEP ply{2}left ply{2}left tra {0} 0 0 {1}\n".format(Lx,divX,n[1]))
        file.write("SETC ply{0}left\n".format(n[1]))   
        file.write("SETR ply{0}left p all\nSETR ply{0}left l all\nSETR ply{0}left s all\n".format(n[1]))     
        for i in range (n[1]-1):
            file.write("SETO ply{0}left\n".format(n[1]-1-i))
            file.write("SWEP plan{0}left plan{1}left tra 0 {2} 0 {3}\n".format(n[1]-i,n[1]-1-i,round(-self.Left.Layup[n[1]-2-i].GetThickness()/W,16),divPly))
            file.write("SWEP ply{0}left ply{0}left tra {1} 0 0 {2}\n".format(n[1]-1-i,Lx,divX))
            file.write("SETC ply{0}left\n".format(n[1]-1-i))   
        # creation of the right down reference point
        file.write("SETO down\nSETO down0\nSETO right\n")
        file.write("PNT p2 {0} {1} {2}\n".format(-0.5*Lx,-0.5+self.OffsetY/W,-0.5*Lz+self.Down.GetTotThickness()/W+self.OffsetZ/W))
        file.write("SETC right\n")
        #creation of the down side
        for i in range(n[3]): #right side
            file.write("SWEP down{2} down{3} tra 0 {0} 0 {1}\n".format(round(self.Right.Layup[n[3]-1-i].GetThickness()/W,16),divPly,i,i+1))
        file.write("SWEP down{2} down{3} tra 0 {0} 0 {1}\n".format(round(1-(self.Left.GetTotThickness()+self.Right.GetTotThickness())/W,16),divYZ,n[3],n[3]+1))
        for i in range(n[1]): #left side
            file.write("SWEP down{2} down{3} tra 0 {0} 0 {1}\n".format(round(self.Left.Layup[i].GetThickness()/W,16),divPly,n[3]+1+i,n[3]+2+i))

        # ~ file.write("SWEP down down tra 0 1 0 {0}\n".format(divYZ))
        file.write("SETO ply1down\n")
        file.write("SWEP down plan1down tra 0 0 {0} {1}\n".format(-self.Down.Layup[0].GetThickness()/W,divPly))
        file.write("SWEP ply1down ply1down tra {0} 0 0 {1}\n".format(Lx,divX))
        file.write("SETR ply1down p all\nSETR ply1down l all\nSETR ply1down s all\n")
        file.write("SETC ply1down\n")        
        for i in range (n[2]-1):
            file.write("SETO ply{0}down\n".format(i+2))
            file.write("SWEP plan{0}down plan{1}down tra 0 0 {2} {3}\n".format(i+1,i+2,round(-self.Down.Layup[i+1].GetThickness()/W,16),divPly))
            file.write("SWEP ply{0}down ply{1}down tra {2} 0 0 {3}\n".format(i+2,i+2,Lx,divX))
            file.write("SETR ply{0}down p all\nSETR ply{0}down l all\nSETR ply{0}down s all\n".format(i+2))
            file.write("SETC ply{0}down\n".format(i+2))        
        # ~ #creation of the right side
        file.write("SETC down\nSETO right\n") 
        file.write("SWEP right right tra 0 0 {0} {1}\n".format(round((H-self.Up.GetTotThickness()-self.Down.GetTotThickness())/W,16),divYZ))
        file.write("SETO ply{0}right\n".format(n[3]))
        file.write("SWEP right plan{2}right tra 0 {0} 0 {1}\n".format(round(self.Right.Layup[n[3]-1].GetThickness()/W,16),divPly,n[3]))
        file.write("SWEP ply{2}right ply{2}right tra {0} 0 0 {1}\n".format(Lx,divX,n[3]))
        file.write("SETC ply{0}right\n".format(n[3]))   
        file.write("SETR ply{0}right p all\nSETR ply{0}right l all\nSETR ply{0}right s all\n".format(n[3]))     
        for i in range (n[3]-1):
            file.write("SETO ply{0}right\n".format(n[3]-1-i))
            file.write("SWEP plan{0}right plan{1}right tra 0 {2} 0 {3}\n".format(n[3]-i,n[3]-1-i,round(self.Right.Layup[n[3]-2-i].GetThickness()/W,16),divPly))
            file.write("SWEP ply{0}right ply{0}right tra {1} 0 0 {2}\n".format(n[3]-1-i,Lx,divX))
            file.write("SETC ply{0}right\n".format(n[3]-1-i))   
        # Merging and meshing
        file.write('MERG p all\n')    
        file.write('MERG l all\n')    
        file.write('MERG s all\n') 
        file.write("ELTY all "+TypeElem+"\n")
        file.write("MESH all\n")
        file.write("SEND all abq\n")
        for i in range (n[0]):
            file.write("SEND ply{0}up abq nam\n".format(i+1))   
        for i in range (n[1]):
            file.write("SEND ply{0}left abq nam\n".format(i+1))   
        for i in range (n[2]):
            file.write("SEND ply{0}down abq nam\n".format(i+1))   
        for i in range (n[0]):
            file.write("SEND ply{0}right abq nam\n".format(i+1))   
        file.close
        
        
    ## Create the input file for ccx solver
    #@param Stress True = output the stress tensor for each elementary load case, False = no stress output
    #@param PlaneSection True = warping of the cross section is not allowed, False = warping of the cross section is allowed
    #@param Disp  Determine the set of elementary load cases : 0=all the 4 load cases, 1 = traction, 2 = warping, 3 = bending span-wise, 4 = bending chord-wise.
    def CreateInpFile(self,Stress=False,PlaneSection=False,Disp=0):
        # 0 = Up; 1 = Left; 2 = Down; 3 = Right    
        Materials = self.Up.Materials
        Orientations = self.Up.Orientations
        NbPly = []
        NbPly.append(len(self.Up.Layup))
        NbPly.append(len(self.Left.Layup))
        NbPly.append(len(self.Down.Layup))
        NbPly.append(len(self.Right.Layup))
        for i in range(len(self.Left.Materials)):
            if self.Left.Materials[i] not in Materials:
                Materials.append(self.Left.Materials[i])
        for i in range(len(self.Left.Orientations)):
            if self.Left.Orientations[i] not in Orientations:
                Orientations.append(self.Left.Orientations[i]) 
        for i in range(len(self.Down.Materials)):
            if self.Down.Materials[i] not in Materials:
                Materials.append(self.Down.Materials[i])
        for i in range(len(self.Down.Orientations)):
            if self.Down.Orientations[i] not in Orientations:
                Orientations.append(self.Down.Orientations[i])  
        for i in range(len(self.Right.Materials)):
            if self.Right.Materials[i] not in Materials:
                Materials.append(self.Right.Materials[i])
        for i in range(len(self.Right.Orientations)):                
            if self.Right.Orientations[i] not in Orientations:
                Orientations.append(self.Right.Orientations[i])  
        file=open("Box.inp","w")
        file.write("*include,input=all.msh\n")
        for i in range(NbPly[0]):
            file.write("*include,input=ply{0}up.nam\n".format(i+1))
        for i in range(NbPly[1]):
            file.write("*include,input=ply{0}left.nam\n".format(i+1))
        for i in range(NbPly[2]):
            file.write("*include,input=ply{0}down.nam\n".format(i+1))
        for i in range(NbPly[3]):
            file.write("*include,input=ply{0}right.nam\n".format(i+1))
        file.write("*include,input=periodic.equ\n")    
        # materials
        i=0
        for Mat in Materials:
            i=i+1
            if isinstance(Mat,OrthoMaterial):
                El = Mat.GetOrtho()[0]
                Et = Mat.GetOrtho()[1]
                Nult = Mat.GetOrtho()[2]
                Nutl = round(float(Et)/El*Nult,6)
                Glt = Mat.GetOrtho()[3]
                Gtl = round(float(Et)/(2+2*Nutl),0)
                file.write("*material,name=MAT{0}\n".format(i))
                file.write("*elastic,type=ENGINEERING CONSTANTS\n")
                file.write("{0},{1},{2},{3},{4},{5},{6},{7}\n{8}\n".format(El,Et,Et,Nult,Nult,Nutl,Glt,Glt,Gtl))
            elif isinstance(Mat,IsoMaterial):
                E = Mat.GetIso()[0]
                Nu = Mat.GetIso()[1]
                file.write("*material,name=MAT{0}\n".format(i))
                file.write("*elastic,type=ISO\n")
                file.write("{0},{1}\n".format(E,Nu))
        # ~ # orientation
        i=0
        for Ori in Orientations:
            i=i+1
            alpha = np.pi/180*Ori
            cos = round(math.cos(alpha),6)
            sin = round(math.sin(alpha),6)
            file.write("*orientation,name=Ori{0}up\n".format(i))
            file.write("{0},{1},{2},{3},{4},{5}\n".format(cos,sin,0.,-sin,cos,0.))
            file.write("*orientation,name=Ori{0}left\n".format(i))
            file.write("{0},{1},{2},{3},{4},{5}\n".format(cos,0.,sin,-sin,0.,cos))
            file.write("*orientation,name=Ori{0}down\n".format(i))
            file.write("{0},{1},{2},{3},{4},{5}\n".format(cos,-sin,0.,-sin,-cos,0.))
            file.write("*orientation,name=Ori{0}right\n".format(i))
            file.write("{0},{1},{2},{3},{4},{5}\n".format(cos,0.,-sin,-sin,0.,-cos))
        # elements affectation
        for i in range(NbPly[0]): #side up
            iMat = Materials.index(self.Up.Layup[i].GetMaterial())+1
            iOri = Orientations.index(self.Up.Layup[i].GetOrientation())+1
            file.write("*solid section, material=MAT{0},elset=Eply{1}up".format(iMat,i+1))
            if isinstance(Materials[iMat-1],OrthoMaterial):
                file.write(",orientation=Ori{0}up".format(iOri))
            file.write("\n")
        for i in range(NbPly[1]): #side left
            iMat = Materials.index(self.Left.Layup[i].GetMaterial())+1
            iOri = Orientations.index(self.Left.Layup[i].GetOrientation())+1
            file.write("*solid section, material=MAT{0},elset=Eply{1}left".format(iMat,i+1))
            if isinstance(Materials[iMat-1],OrthoMaterial):
                file.write(",orientation=Ori{0}left".format(iOri))
            file.write("\n")
        for i in range(NbPly[2]): #side down
            iMat = Materials.index(self.Down.Layup[i].GetMaterial())+1
            iOri = Orientations.index(self.Down.Layup[i].GetOrientation())+1
            file.write("*solid section, material=MAT{0},elset=Eply{1}down".format(iMat,i+1))
            if isinstance(Materials[iMat-1],OrthoMaterial):
                file.write(",orientation=Ori{0}down".format(iOri))
            file.write("\n")
        for i in range(NbPly[3]): #side right
            iMat = Materials.index(self.Right.Layup[i].GetMaterial())+1
            iOri = Orientations.index(self.Right.Layup[i].GetOrientation())+1
            file.write("*solid section, material=MAT{0},elset=Eply{1}right".format(iMat,i+1))
            if isinstance(Materials[iMat-1],OrthoMaterial):
                file.write(",orientation=Ori{0}right".format(iOri))
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
        
        
    ## Compute the CrossSection::MassMatrix using the mass matrix of the 4 walls of the box
    #@return Mu Mass per unit length (kg/m)
    #@return I22 mass moment of inertia around yB (bending span-wise, kg.m) 
    #@return I33 mass moment of inertia around zB (bending chord-wise, kg.m) 
    #@return I23 product of inertia in plane (yB,zB)
    #@return Ycg Y coordinate of the Center of Gravity in frame B
    #@return Zcg Z coordinate of the Center of Gravity in frame B
    def ComputeMassMatrix(self,OffsetY=0.,OffsetZ=0.):
        OffsetYup = OffsetY
        OffsetZup = OffsetZ - 0.5*self.Height + 0.5*self.Up.GetTotThickness()
        OffsetYleft = -OffsetZ
        OffsetZleft = OffsetY - 0.5*self.Width + 0.5*self.Left.GetTotThickness()
        OffsetYdown = OffsetY
        OffsetZdown = OffsetZ + 0.5*self.Height - 0.5*self.Down.GetTotThickness()
        OffsetYright = -OffsetZ
        OffsetZright = OffsetY  + 0.5*self.Width - 0.5*self.Right.GetTotThickness()
        #side up
        [Muup,I22up,I33up,I23up,Ycgup,Zcgup]=self.Up.ComputeMassMatrix(OffsetY=OffsetYup,OffsetZ=OffsetZup)
        [Muleft,I33left,I22left,I23left,Zcgleft,Ycgleft]=self.Left.ComputeMassMatrix(OffsetY=OffsetYleft,OffsetZ=OffsetZleft)
        [Mudown,I22down,I33down,I23down,Ycgdown,Zcgdown]=self.Down.ComputeMassMatrix(OffsetY=OffsetYdown,OffsetZ=OffsetZdown)
        [Muright,I33right,I22right,I23right,Zcgright,Ycgright]=self.Right.ComputeMassMatrix(OffsetY=OffsetYright,OffsetZ=OffsetZright)
        Zcgleft = -Zcgleft
        Zcgright = -Zcgright
        I23left = - I23left
        I23right = -I23right
        Mu = Muup+Muleft+Mudown+Muright
        Ycg = (Muup*Ycgup+Muleft*Ycgleft+Mudown*Ycgdown+Muright*Ycgright)/Mu
        Ycg = round(1+Ycg,12)-1        
        Zcg = (Muup*Zcgup+Muleft*Zcgleft+Mudown*Zcgdown+Muright*Zcgright)/Mu
        Zcg = round(1+Zcg,12)-1
        I22 = I22up+I22left+I22down+I22right
        I33 = I33up+I33left+I33down+I33right
        I23 = I23up+I23left+I23down+I23right
        return [Mu,I22,I33,I23,Ycg,Zcg]
        
        
    ## Interface to the utils subroutine utils::CreatePeriodicEq
    def CreatePeriodicEq(self):
        [self.nstrain,self.ncurv,self.x0,self.Lx,self.nodes,self.elements]=CreatePeriodicEq("all.msh")
    
    ## Launch the cgx postprocessor with a particular elementary load case.
    #@param TypeElem the type of finite element used for the computation. Supported : he20r, he20, he8i, he8, pe15
    #@param NbElemX the number of finite element in the beam direction (for a constant cross section, 1 element is enough)
    #@param NBElemYZ the number of finite element along the wall directions
    #@param NBElemPly the number of finite element in a composite ply
    #@param DefType  Determine the set of elementary load cases : 0=all the 4 load cases, 1 = traction, 2 = warping, 3 = bending span-wise, 4 = bending chord-wise.
    #@param PlaneSection : True = warping of the cross section is not allowed, False = warping of the cross section is allowed
    def DisplaySectionDeformation(self,TypeElem,NbElemX,NbElemYZ,NbElemPly,DefType,PlaneSection=False):
        self.CreateFbdFile(TypeElem,NbElemX,NbElemYZ,NbElemPly)
        RunFbdFile("Box.fbd")
        self.CreatePeriodicEq()
        self.CreateInpFile(Disp=DefType,PlaneSection=PlaneSection)
        Status = RunInpFile("Box")[0]
        if Status !=0:
            raise RuntimeError("calculix computation failure")
        RunFrdFile("Box.frd")
