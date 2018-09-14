# coding=UTF-8

## @package utils The package contains a set of functions used by several classes.
import numpy as np
import subprocess as sp
import math
import cmath
                
## This function create the linear equations between nodes ddl from opposite periodic faces using the ccx input card *EQUATION (see GEBTAero paper for more details).
#@param Meshfile the name of the MeshFile used in the analysis (either creating using cgx or external mesh)
#@param OffsetY Y coordinate of the box center in the Frame B
#@param OffsetZ Z coordinate of the box center in the Frame B
#@return nstrain index of the nstrain dummy node
#@return ncurv index of the ncurv dummy node
#@return x0 list of node index belonging to the reference periodic RVE face
#@return Lx Length of the RVE in the beam axis direction
#@return na list of nodes (INP format)
#@return elements list of elements (INP format)
def CreatePeriodicEq(MeshFile,OffsetY=0.,OffsetZ=0.):
    # initializing
    datatype="unknown"
    nodes=[]
    elements=[]
    IndElt = []
    elsets=[]
    f = open(MeshFile,"r")
    fo = open("periodic.equ","w")

    # read node table
    for line in f:
        if line.startswith("*"):
            # remove spaces and make uppercase
            line=(line.upper()).replace(" ","")
            if line.startswith("*NODE"):
                datatype="node"
            elif line.startswith("*ELEMENT"):  
                datatype="element"               
            elif line.startswith("*ELSET"):  
                datatype="elset"
                temp = line
                temp = temp.split('=')
                ElsetName = temp[1].replace("\n","")            
            else:
                datatype="unknown"
            continue    
        if datatype=="node":
            temp = [float(field) for field in line.split(",")]
            temp[0] = int(temp[0])
            nodes.append(temp)    
        elif datatype=="element":
            temp = []
            line = line.split(",")
            if len(line) >2:
                IndElt.append(int(line[0]))
                for i in range(len(line)-1):
                    temp.append(int(line[i]))
                elements.append(temp)
        elif datatype == "elset":
            line = line.split(",")
            for i in range(len(line)-1):
                elt = int(line[i])
                ind = IndElt.index(elt)
                elements[ind].append(ElsetName)                
   
    f.close()
    na=np.array(nodes)
    n=np.int_(na[:,0])
    x=np.array(na[:,1])
    y=np.array(na[:,2])
    z=np.array(na[:,3])

    # determine dimensions of the brick
    [lx0,ly0,lz0]=[x.min(),y.min(),z.min()]
    [lx1,ly1,lz1]=[x.max(),y.max(),z.max()]
    [Lx,Ly,Lz] = [lx1-lx0,ly1-ly0,lz1-lz0]
    xB = 0.5*(lx0+lx1)
    # determine control nodes
    x0 = []
    xl = []
    for i in range(len(na)):
        if math.isclose(x[i],lx0,abs_tol=1e-12):
            x0.append(n[i])
        elif math.isclose(x[i],lx1,abs_tol=1e-12):
            xl.append(n[i])         
          
    """
    Create the equations:
    ux_j=ux_i+ux_nx on all inner points
    """
    # create a dummy node containing the strain curvature 
    ncurv = max(n)+1 
    nstrain = ncurv+1       
    fo.write("*node,Nset=ncurv\n{0},{1},{2},{3}\n".format(ncurv,round(0.5*(lx0+lx1),8),round(0.5*(ly0+ly1),8),+round(0.5*(lz0+lz1),8)))  
    fo.write("*node,Nset=nstrain\n{0},{1},{2},{3}\n".format(nstrain,round(0.5*(lx0+lx1),8),round(0.5*(ly0+ly1),8),round(0.5*(lz0+lz1),8)))  
    
    # print x1 node set
    fo.write("*Nset,NSET=nxl\n")
    for node in xl:
        fo.write("{0},\n".format(node))
    fo.write("\n")           
        
    # print x0 node set
    fo.write("*Nset,NSET=nx0\n")
    for node in x0:
        fo.write("{0},\n".format(node))
    fo.write("\n")           
        
    fo.write("*equation\n")
   
    # eliminate nodes at x0
    for i in x0:
        for j in xl:
            if (y[i-1]==y[j-1] and z[i-1]==z[j-1]):
                #~ y1=y[j-1]-0.5*(ly0+ly1)-OffsetY
                #~ z1=z[j-1]-0.5*(lz0+lz1)-OffsetZ
                y1=y[j-1]-OffsetY
                z1=z[j-1]-OffsetZ
                # -ux0 + ux1 - Lx*Epsi - Lx*z1*Gamma2 + Lx*y1*Gamma3
                fo.write("5\n{0},1,-1,{1},1,1,{2},1,{3},{4},2,{5},{6},3,{7}\n".format(i,j,nstrain,round(-Lx,6),ncurv,round(-Lx*z1,6),ncurv,round(Lx*y1,6)))
                # -uy0 + uy1 + xB*Lx*Gamma3 + Lx*z1*Gamma1
                fo.write("4\n{0},2,-1,{1},2,1,{2},3,{3},{4},1,{5}\n".format(i,j,ncurv,round(xB*Lx,6),ncurv,round(Lx*z1,6)))
                # -uz0 + uz1 + xB*Lx*Gamma2 -Lx*y1*Gamma1
                fo.write("4\n{0},3,-1,{1},3,1,{2},2,{3},{4},1,{5}\n".format(i,j,ncurv,round(xB*Lx,6),ncurv,round(-Lx*y1,6)))
                # ~ # -ux0 + ux1 - Lx*Epsi - Lx*z1*Gamma2 + Lx*y1*Gamma3
                # ~ fo.write("5\n{0},1,-1,{1},1,1,{2},1,{3},{4},2,{5},{6},3,{7}\n".format(i,j,nstrain,1.,ncurv,round(-z1,6),ncurv,round(y1,6)))
                # ~ # -uy0 + uy1 + xB*Lx*Gamma3 + Lx*z1*Gamma1
                # ~ fo.write("4\n{0},2,-1,{1},2,1,{2},3,{3},{4},1,{5}\n".format(i,j,ncurv,round(xB,6),ncurv,round(z1,6)))
                # ~ # -uz0 + uz1 + xB*Lx*Gamma2 -Lx*y1*Gamma1
                # ~ fo.write("4\n{0},3,-1,{1},3,1,{2},2,{3},{4},1,{5}\n".format(i,j,ncurv,round(xB,6),ncurv,round(-y1,6)))
    return [nstrain,ncurv,x0,Lx,na,elements]

## launch cgx preprocessor to create a .msh mesh file
#@param FileName The name of the input file
def RunFbdFile(FileName):
    Command = "cgx -bg "+FileName
    return sp.getstatusoutput(Command)

## Launch ccx solver to create a .dat file
#@param FileName The name of the input file
def RunInpFile(FileName):
    Command = "ccx "+FileName
    return sp.getstatusoutput(Command)

## Launch cgx postprossor
#@param FileName The name of the input file
def RunFrdFile(FileName):
    Command = "cgx "+FileName
    return sp.getstatusoutput(Command)
    
## Launch a Paraview script    
#@param FileName The name of the input file
def RunParaviewScript(FileName):
    Command = "paraview --script="+FileName    
    return sp.getstatusoutput(Command)
    
## Launch the unical mesh converter (from unv to inp format)  
#@param FileName The name of the input file
#@param FileOut the name of the ouptut mesh file
def RunUnvConv(FileName,FileOut):
    Command = "unical "+FileName+" "+FileOut
    return sp.getstatusoutput(Command)   
    
## Remove of file generated by gebtaero, ccx and cgx. 
def RemoveFiles():
    Command = "rm *.nam *.fbd *.inp *.out *.dat *.equ *.cvg *.sta *.frd *.ech *.ini *.vec"
    return sp.getstatusoutput(Command)
    
## Romeve .msh files
def RemoveMeshFiles():
    Command = "rm *.msh"
    return sp.getstatusoutput(Command)
    
##  Read node field in a .frd file generated by ccx
#@param FileName The name of the input file
#@param FieldName The name of the field as it is written in the .frd file
#@return FieldList A list containing the stored fielddata
def ReadNodesField(self,FileName,FieldName):
    file = open(FileName,"r")
    lines = file.readlines()
    FieldIndex = []
    FieldList = []
    
    for i in range(len(lines)):
        if lines[i].startswith(" -4"):
            FieldIndex.append(i)
    for i in FieldIndex:
        Field = []
        if lines[i].split()[1] == FieldName:
            j = i+1
            while lines[j].startswith(" -5"):
                j = j+1
            while lines[j].startswith(" -1"):
                n = len(lines[j])
                FieldLine = []
                if n < 51:
                    FieldLine.append(int(lines[j][4:13]))
                    FieldLine.append(float(lines[j][13:25]))
                    FieldLine.append(float(lines[j][25:37]))
                    FieldLine.append(float(lines[j][37:49]))
                else:
                    FieldLine.append(int(lines[j][4:13]))
                    FieldLine.append(float(lines[j][13:25]))
                    FieldLine.append(float(lines[j][25:37]))
                    FieldLine.append(float(lines[j][37:49]))
                    FieldLine.append(float(lines[j][49:61]))
                    FieldLine.append(float(lines[j][61:73]))
                    FieldLine.append(float(lines[j][73:85]))
                Field.append(FieldLine)
                j = j+1
            FieldList.append(Field)
    return FieldList
    
## This routine read the .dat file generated by ccx solver and generate the associated CrossSection::FlexibilityMatrix
#@param FileName The name of the input file
#@param nstrain The index of the nstrain dummy node    
#@param ncurv The index of the ncurv dummy node   
#@param Lx Length of the RVE in the beam direction 
#@param tol set to zero value < tol
#@param RigidX if True suppress the traction ddl
#@param RigidZ if True suppress the bending in the lag axis ddl
#@return the CrossSection::FlexibilityMatrix
def ReadFlexibilityFromDisp(FileName,nstrain,ncurv,Lx,tol,RigidX=False,RigidZ=False):
    FlexMat = np.zeros([6,6])
    StrainIndex = []
    CurvIndex = []
    file = open(FileName,"r")
    lines = file.readlines()
    n = len(lines)
    for i in range(n):
        temp = lines[i].split()
        if len(temp) > 0:
            if temp[0] == str(nstrain):
                StrainIndex.append(i)
            if temp[0] == str(ncurv):
                CurvIndex.append(i)          
    if len(StrainIndex) !=4 or len(CurvIndex) !=4:
        raise RuntimeError("Unable to read all the flexibility matrix coefficients")      
    # first column (traction)
    FlexMat[0,0] = Lx*float(lines[StrainIndex[0]].split()[1])
    FlexMat[3,0] = FlexMat[0,3] = 0.5*Lx*(float(lines[StrainIndex[1]].split()[1])+float(lines[CurvIndex[0]].split()[1]))
    FlexMat[4,0] = FlexMat[0,4] = 0.5*Lx*(float(lines[StrainIndex[2]].split()[1])+float(lines[CurvIndex[0]].split()[2]))
    FlexMat[5,0] = FlexMat[0,5] = 0.5*Lx*(float(lines[StrainIndex[3]].split()[1])+float(lines[CurvIndex[0]].split()[3]))
    # column four (twisting)
    FlexMat[3,3] = Lx*float(lines[CurvIndex[1]].split()[1])        
    FlexMat[4,3] = FlexMat[3,4] = 0.5*Lx*(float(lines[CurvIndex[2]].split()[1])+float(lines[CurvIndex[1]].split()[2]))      
    FlexMat[5,3] = FlexMat[3,5] = 0.5*Lx*(float(lines[CurvIndex[3]].split()[1])+float(lines[CurvIndex[1]].split()[3]))     
    # column five (bending spanwise)
    FlexMat[4,4] = Lx*float(lines[CurvIndex[2]].split()[2])        
    FlexMat[5,4] = FlexMat[4,5] = 0.5*Lx*(float(lines[CurvIndex[3]].split()[2])+float(lines[CurvIndex[2]].split()[3]))           
    # column six (bending chordwise)
    FlexMat[5,5] = Lx*float(lines[CurvIndex[3]].split()[3])   
    
    # suppression of first column/line if RigidX = True (no axial strain)
    if (RigidX):
        FlexMat[0,:] = FlexMat[:,0] = 0.
    # suppression of last column/line if RigidZ=True (no lag ddl)
    if (RigidZ):
        FlexMat[5,:] = FlexMat[:,5] = 0.
    # suppression of numerical noise
    MaxFlex = np.max(abs(FlexMat))
    for i in range(6):
        for j in range(6):
            if abs(FlexMat[i,j]) < tol*MaxFlex:
                FlexMat[i,j] = 0.         
    return FlexMat
    
## This function compute the correlation coefficients between two matrices (used to follow an aeroelastic mode between successive upstream flow velocity)
#@param Tab1 Matrix 1
#@param Tab2 Matrix 2
#@param Index the correlation computation is restricted to the line number in Index (to the followed aeroelastic mdoes)
#@return Corr The correlation matrix between Tab1 and Tab2
def CorrelateTab(Tab1,Tab2,Index):
    n = len(Tab1)
    m = len(Tab2)
    Corr = np.zeros([n,m,2])
    for i in range(n):
        if i in Index:
            for j in range(m):
                #Module correlation
                veca = np.array(Tab1[i][1])
                veca = veca/np.linalg.norm(veca)
                vecb = np.array(Tab2[j][1])
                vecb = vecb/np.linalg.norm(vecb)
                Corr[i,j,0] = abs(np.correlate(veca,vecb))
                # Eigenvalue correlation
                veca = np.array(Tab1[i][0])
                veca = veca/np.linalg.norm(veca)
                vecb = np.array(Tab2[j][0])
                vecb = vecb/np.linalg.norm(vecb)
                Corr[i,j,1] = abs(np.correlate(veca,vecb))
                #~ #Phase correlation
                #~ veca = np.array(Tab1[i][1])*np.array(Tab1[i][2])
                #~ veca = veca/np.linalg.norm(veca)
                #~ vecb = np.array(Tab2[j][1])*np.array(Tab2[j][2])
                #~ vecb = vecb/np.linalg.norm(vecb)
                #~ Corr[i,j] = Corr[i,j]*abs(np.correlate(veca,vecb))
    return Corr
    
## Compute the position of the Elastic Center using CrossSection::FlexibilityMatrix coeffcients
#@param FlexMat the input  CrossSection::FlexibilityMatrix
#@return ECy the Y coordinate of the Elastic Center in Frame B
#@return ECz the Z coordinate of the Elastic Center in Frame B
def ComputeElasticCenterFromFlexMat(FlexMat):
    if(math.isclose(FlexMat[4,4],0.) or math.isclose(FlexMat[5,5],0.)):
        raise RuntimeError("the curvature flexibility must be nonzero to compute the elastic center position")
    else :
        ECy = float(FlexMat[5,0])/FlexMat[5,5]
        ECz = -float(FlexMat[4,0])/FlexMat[4,4]
    return [ECy,ECz]
    
## Launch gebtaero using Command and read the output values in the shell  
#@param Command the shell command to launch gebtaero 
#@return  EigenValues A list containing the eigenvalues
#@return  EigenVectors A list containing the eigenvectors
#@return  StaticLoads A list containing wingroot static loads
def ReadFromPipe(Command):
    EigenValues = []
    EigenVectors = []
    StaticLoads = []
    Vector = []
    result= sp.Popen(Command,stdout=sp.PIPE,stderr=sp.PIPE)
    out, err = result.communicate()
    errcode = result.returncode
    #~ # search for gebtaero error
    if errcode !=0:
        print(err, end="\n") 
        raise RuntimeError("error in gebtaero solver") 
    lines = [str(field).replace("b'","") for field in out.splitlines()]
    datatype="unknown"
    for line in lines:
        if line.strip().startswith("*"):
            # remove spaces and make uppercase
            line=(line.upper()).replace(" ","")
            if line.startswith("*EIGENVALUES"):
                datatype="eigenvalues"
            elif line.startswith("*EIGENVECTORS"):  
                datatype="eigenvectors"                
            elif line.startswith("*STATICLOADS"):  
                datatype="staticloads"                
            else:
                datatype="unknown"
            continue    
        if datatype=="eigenvalues":
            temp = [float(field) for field in line.replace("'","").split()]
            EigenValues.append(temp)    
        elif datatype=="eigenvectors":
                temp = [float(field) for field in line.replace("'","").split()]
                EigenVectors.append(temp)
        elif datatype=="staticloads":
            temp = [float(field) for field in line.replace("'","").split()]
            StaticLoads.append(temp)
    return EigenValues,EigenVectors,StaticLoads  
  
## Launch gebtaero aero in the shell with Command and get back eigenproblem data
#@param Command the shell command to launch gebtaero 
#@param output The type of data to output. output = modes : eigenvalues and eigenvectors, output = eigenvalues : eigenvalues only.
#@return EigenModes a list of modes containing for each of them the eigenvalue and the eigenvector
#@return EigenValues a list of eigenvalues
def ReadModesFromPipe(Command,output="modes"):
    EigenModes = []
    mode = []
    counter = 0
    result= sp.Popen(Command,stdout=sp.PIPE,stderr=sp.PIPE)
    out, err = result.communicate()
    errcode = result.returncode
    #~ # search for gebtaero error
    if errcode == 106:
        pass
    elif errcode !=0:
        print(err, end="\n") 
        raise RuntimeError("error in gebtaero solver") 
    lines = [str(field).replace("b'","") for field in out.splitlines()]
    datatype="unknown"
    for line in lines:
        if line.strip().startswith("*"):
            # remove spaces and make uppercase
            line=(line.upper()).replace(" ","")
            if line.startswith("*REAL"):
                datatype="real"
            elif line.startswith("*COMPLEX"):  
                datatype="complex"                
            else:
                datatype="unknown"
            continue    
        if datatype=="real":
            mode.append([float(field) for field in line.replace("'","").split()])
            counter = counter+1
            if output == "modes":
                if counter >=2:
                    mode.append([0.]*len(mode[1]))
                    for j in range(len(mode[1])):
                        vec = []
                        temp = complex(mode[1][j],0.)
                        vec.append(np.abs(temp))
                        vec.append(np.angle(temp))
                        mode[1][j]=vec[0]
                        mode[2][j]=vec[1]
                    EigenModes.append(mode)
                    mode = []
                    counter = 0
            else:
                if counter >=1:            
                    EigenModes.append(mode)
                    mode = []
                    counter = 0
        elif datatype=="complex":
            mode.append([float(field) for field in line.replace("'","").split()])
            counter = counter+1
            if output == "modes":
                if counter >=3:
                    for j in range(len(mode[1])):
                        vec = []
                        temp = complex(mode[1][j],mode[2][j])
                        vec.append(np.abs(temp))
                        vec.append(np.angle(temp))
                        mode[1][j]=vec[0]
                        mode[2][j]=vec[1]
                    # ~ mode[1][0] = math.sqrt(mode[1][0]*mode[2][0])
                    EigenModes.append(mode)
                    mode = []
                    counter = 0
            else:
                if counter >=1:    
                    EigenModes.append(mode)
                    mode = []
                    counter = 0    

    if output == "modes":       
        return EigenModes,errcode
    elif output == "eigenvalues":  
        EigenValues = []
        for i in range(len(EigenModes)):
            EigenValues.append(EigenModes[i][0])
        return EigenValues    
  
## Launch gebtaero in the shell with Command and get back the wing root static loads
#@param Command the shell command to launch gebtaero 
#@param output The type of data to output. output="static" : the wing root static loads
#@return Loads a list 2*3 containing wing root forces and moments
def ReadLoadsFromPipe(Command,output="static"):
    Loads = []
    result= sp.Popen(Command,stdout=sp.PIPE,stderr=sp.PIPE)
    out, err = result.communicate()
    errcode = result.returncode
    #~ # search for gebtaero error
    if errcode !=0:
        print(err, end="\n") 
        raise RuntimeError("error in gebtaero solver") 
    lines = [str(field).replace("b'","") for field in out.splitlines()]
    datatype="unknown"
    for line in lines:
        if line.strip().startswith("*"):
            # remove spaces and make uppercase
            line=(line.upper()).replace(" ","")
            if line.startswith("*STATIC"):
                datatype="static"             
            else:
                datatype="unknown"
            continue    
        if datatype=="static":
            temp = [float(field) for field in line.replace("'","").split()]
            Loads.append(temp)
    if output == "static":       
        return Loads
  
