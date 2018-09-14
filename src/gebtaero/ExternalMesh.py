# coding=UTF-8
import numpy as np
import math
import os
from .utils import *
from .OrthoMaterial import *
from .IsoMaterial import *

##This class allow to use a mesh created without cgx, the elset name must correspond to the component name. The input mesh format has to be inp format or unv using the unical format translator
class ExternalMesh:
    
    def __init__(self,MeshFile,OffsetY=0.,OffsetZ=0.,UnvConv=False,Chord=1.):
        # does the meshfile exist
        if not os.path.isfile(MeshFile):
            raise FileNotFounfError("The mesh file is not found")
        if UnvConv:            
            RunUnvConv(MeshFile,"converted.msh")
            self.MeshFile = "converted.msh"
        else:
            ## the name of the meshfile
            self.MeshFile = MeshFile    
        ## a list of group of elements (material affectation)
        self.Components = []
        ## a list of materials (either IsoMaterial or OrthoMaterial)
        self.Materials = []
        ## a list of fiber orientation
        self.Orientations = []
        ## Y coordinate of the mesh center in the Frame B
        self.OffsetY = OffsetY
        ## Z coordinate of the mesh center in the Frame B
        self.OffsetZ = OffsetZ
        ## The chord of the airfoil
        self.Chord=Chord
        [nstrain,ncurv,x0,Lx,nodes,elements]=CreatePeriodicEq(self.MeshFile)
        RemoveFiles()
        ## index of the nstrain ccx dummy node
        self.nstrain = None
        ## index of the nstrain ccx dummy node
        self.ncurv = None
        ## length of the RVE in direction xb
        self.Lx = None
        ## list of the mesh nodes (inp format)
        self.nodes = None
        ## index of the nodes belonging to the first RVE periodic face 
        self.x0 = None
        ## list of the mesh element (inp format)
        self.elements = None
        
    ## interface to the utils::CreatePeriodicEq
    def CreatePeriodicEq(self):
        [self.nstrain,self.ncurv,self.x0,self.Lx,self.nodes,self.elements]=CreatePeriodicEq(self.MeshFile)            
        
    def GetMeshFile(self):
        return self.MeshFile    
        
    ## Link a mesh element group to a material
    #@param Name the name of the element group as it is set in the meshfile (ELSET)
    #@param Material material to be set in the element group (either IsoMatrial or OrthoMaterial)
    #@param the fiber orientation in case of OrthoMaterial
    def AppendComponent(self,Name,Material,Orientation=None):
        self.Components.append([Name.upper(),Material,Orientation])
        if(Material not in self.Materials):
            self.Materials.append(Material)
        if Orientation is None :
            if isinstance(Material,OrthoMaterial):
                raise RuntimeError("An orientation must be defined for a composite material")
        else:    
            if(Orientation not in self.Orientations):
                self.Orientations.append(Orientation)
            
    ## Create the input file for ccx solver
    #@param Stress True = output the stress tensor for each elementary load case, False = no stress output
    #@param PlaneSection True = warping of the cross section is not allowed, False = warping of the cross section is allowed
    #@param Disp  Determine the set of elementary load cases : 0=all the 4 load cases, 1 = traction, 2 = warping, 3 = bending span-wise, 4 = bending chord-wise.        
    def CreateInpFile(self,Stress=False,PlaneSection=False,Disp=0):
        NbComp = len(self.Components)
        NbMat = len(self.Materials)
        NbOri = len(self.Orientations)
        if (NbComp <1):
            raise RuntimeError("You must define at least one component for your mesh")
        file=open("Ext.inp","w")
        file.write("*include,input={0}\n".format(self.MeshFile))
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
                file.write("*material,name=MAT{0}\n".format(i+1))
                file.write("*elastic,type=ENGINEERING CONSTANTS\n")
                file.write("{0},{1},{2},{3},{4},{5},{6},{7}\n{8}\n".format(El,Et,Nult,Nutl,Nutl,Glt,Glt,Gtl))
            elif isinstance(Mat,IsoMaterial):
                E = Mat.GetIso()[0]
                Nu = Mat.GetIso()[1]
                file.write("*material,name=MAT{0}\n".format(i+1))
                file.write("*elastic,type=ISO\n")
                file.write("{0},{1}\n".format(E,Nu))
        # orientation
        for i in range(NbOri):
            alpha = np.pi/180*self.Orientations[i]
            cos = str(round(math.cos(alpha),8))
            sin = str(round(math.sin(alpha),8))
            file.write("*orientation,name=Ori{0}\n".format(i+1))
            file.write("{0},{1},0.,{2},{3},0.\n".format(cos,sin,-sin,cos))
        # elements affectation
        for i in range(NbComp):
            Name = self.Components[i][0]
            iMat = self.Materials.index(self.Components[i][1])+1
            if isinstance(self.Materials[iMat-1],OrthoMaterial):
                iOri = self.Orientations.index(self.Components[i][2])+1
                file.write("*solid section, material=MAT{0},elset={1},orientation=Ori{2}\n".format(iMat,Name,iOri))
            elif isinstance(self.Materials[iMat-1],IsoMaterial):
               file.write("*solid section, material=MAT{0},elset={1}\n".format(iMat,Name))
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
        
    ## Compute the surface and the center of gravity of the triangular face of a penta element
    #@param nodes The node list of the Mesh
    #@param element the index of the element computed
    #@param x0 the nodes index list of the RVE periodic face
    #@return surface the area of the triangular face belonging to the reference periodic face of the RVE (which nodes are in x0)
    #@return ycg the y position of the element CG    
    #@return zcg the z position of the element CG    
    #@return elI22 inertia I22 of the element   
    #@return elI33 inertia I33 of the element   
    #@return elI23 inertia I23 of the element   
    def ComputeElementSurfAndCG(self,nodes,element,x0):
        vertex = []
        # extractions of node belonging to the x0 face
        for i in range(len(element)-2):
            node = element[i+1]
            if  node in x0:
                ind = np.argwhere(nodes==node)[0,0]
                vertex.append(nodes[ind,2:])
        if len(vertex) == 3:
            eltype = 'tri'
        elif len(vertex) == 4:
            eltype = 'quad' 
        else:
            raise RuntimeError("the type of element is not supported for inertia computation.")
        # surface calculation
        if eltype == 'tri':
            ya = vertex[0][0]-self.OffsetY
            za = vertex[0][1]-self.OffsetZ
            yb = vertex[1][0]-self.OffsetY
            zb = vertex[1][1]-self.OffsetZ
            yc = vertex[2][0]-self.OffsetY
            zc = vertex[2][1]-self.OffsetZ
            surface = 0.5*abs((yb-ya)*(zc-za)-(yc-ya)*(zb-za))
            # surface center
            ycg = 1./3*(ya+yb+yc)    
            zcg = 1./3*(za+zb+zc)
            # inertia contribution
            elI22 = ya**2+yb**2+yc**2+ycg**2
            elI33 = za**2+zb**2+zc**2+zcg**2
            elI23 = ya*za+yb*zb+yc*zc+ycg*zcg
        else:
            raise RuntimeError('not coded yet')
        return surface,ycg,zcg,elI22,elI33,elI23

    ## Approximate the CrossSection::MassMatrix using the mesh elements
    #@param nodes The node list of the Mesh
    #@param x0 the nodes index list of the RVE periodic face
    #@param elements list of the mesh elements (inp format)
    #@return Mu Mass per unit length (kg/m)
    #@return I22 mass moment of inertia around yB (bending span-wise, kg.m) 
    #@return I33 mass moment of inertia around zB (bending chord-wise, kg.m) 
    #@return I23 product of inertia in plane (yB,zB)
    #@return Ycg Y coordinate of the Center of Gravity in frame B
    #@return Zcg Z coordinate of the Center of Gravity in frame B
    def ComputeMassMatrixFromMesh(self,nodes,x0,elements):
        Ycg = 0.
        Zcg = 0.
        surface = 0.
        Mu = 0.
        I22 = 0.
        I33 = 0.
        I23 = 0.
        CompArray = np.array(self.Components)
        for elt in elements:
            #determination of the element material
            Comp = elt[-1]
            if not isinstance(Comp,str):
                raise RuntimeError("some elements have no material defined ")
            Ind = np.argwhere(CompArray == Comp)
            Mat = CompArray[Ind[0,0,],1]
            Rho = Mat.GetDensity()
            [elsurface,elycg,elzcg,elI22,elI33,elI23]=self.ComputeElementSurfAndCG(nodes,elt,x0)
            # append elemental properties to the global ones
            surface = surface + elsurface
            Mu= Mu + Rho*elsurface
            Ycg = Ycg + Rho*elsurface*elycg
            Zcg = Zcg + Rho*elsurface*elzcg
            #~ I22 = I22 + Rho*elsurface*elI22
            #~ I33 = I33 + Rho*elsurface*elI33
            #~ I23 = I23 + Rho*elsurface*elI23
            I22 = I22 + Rho*elsurface*elzcg**2
            I33 = I33 + Rho*elsurface*elycg**2
            I23 = I23 + Rho*elsurface*elycg*elzcg
        if math.isclose(Mu,0.):
            raise RuntimeError("The mass per unit length must be non zero")
        Ycg = 1./Mu*Ycg
        Zcg = 1./Mu*Zcg
        #~ I22 = 1./20*I22 - MassSurf*Ycg**2
        #~ I33 = 1./20*I33 - MassSurf*Zcg**2
        #~ print(surface,Ycg,Zcg,I22,I33,I23)
        return [Mu,I22,I33,I23,Ycg,Zcg]
        
    ## Launch the cgx postprocessor with a particular elementary load case.
    #@param DefType  Determine the set of elementary load cases : 0=all the 4 load cases, 1 = traction, 2 = warping, 3 = bending span-wise, 4 = bending chord-wise.
    def DisplaySectionDeformation(self,DefType):
        self.CreatePeriodicEq()
        self.CreateInpFile(Disp=DefType)
        RunInpFile("Ext")
        RunFrdFile("Ext.frd")
