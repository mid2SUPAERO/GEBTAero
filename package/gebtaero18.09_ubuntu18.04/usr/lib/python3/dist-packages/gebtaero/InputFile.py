# coding=UTF-8
import subprocess as sp

class InputFile:
    """
    This class is a mirror of the input data file of the solver
    a method is implemented to write the input file using the class argument
    """        
    def __init__(self,Name,AnalysisFlag,AeroFlag,GravFlag,Niter,Nstep,Nvtk,Nev,ACOmegaa,ACOmegaaTFNumber,ACVa,ACVaTFNumber,Wing,Vinf,Rho,AlphaAC,BetaAC,SimuStart,SimuEnd,Xcg=0.):
        self.Name = Name
        self.FileName = Name+".dat"
        self.AnalysisFlag = AnalysisFlag
        self.AeroFlag = AeroFlag
        self.GravFlag = GravFlag
        self.Niter = Niter
        self.Nstep = Nstep
        self.Nvtk = Nvtk
        self.Nev = Nev
        self.ACOmegaa = ACOmegaa
        self.ACOmegaaTFNumber = ACOmegaaTFNumber
        self.ACVa = ACVa
        self.ACVaTFNumber = ACVaTFNumber
        self.Wing = Wing
        self.Vinf = Vinf
        self.Rho = Rho
        self.AlphaAC = AlphaAC
        self.BetaAC = BetaAC
        self.SimuStart = SimuStart
        self.SimuEnd = SimuEnd
        self.Xcg = Xcg
        
        # Initialisation of time function list
        self.TimeFunctions = []
          
    def GetName(self):
        return self.Name
    
    def GetFileName(self):
        return self.FileName    
        
    def GetAnalysisFlag(self):
        return self.AnalysisFlag    
        
    def AppendTimeFunction(self,TimeFunction):
        self.TimeFunctions.append(TimeFunction)
    
    def GetTimeFunction(self,index):
        return self.TimeFunctions[index]
        
    def WriteInputFile(self):
        """
        This method intend to write the input data file of the solver using the class attributes
        """
        nkp = len(self.Wing.GetKpList())
        nmemb = len(self.Wing.GetWingSections())
        nmate = len(self.Wing.GetCrossSections())
        nframe = len(self.Wing.GetFrames())
        ntimefun = len(self.TimeFunctions)
        file = open(self.FileName,"w")
        # first line of the solver
        file.write(str(self.AnalysisFlag))
        file.write(" "+str(self.AeroFlag))
        file.write(" "+str(self.GravFlag))
        file.write(" "+str(self.Niter))
        file.write(" "+str(self.Nstep))
        file.write(" "+str(self.Nvtk)+"\n")
        # Omegaa and Va
        for i in range(3):
            file.write(str(self.ACOmegaa[i])+" ")
        file.write("\n")
        for i in range(3):
            file.write(str(self.ACOmegaaTFNumber)+" ")
        file.write("\n")
        for i in range(3):
            file.write(str(self.ACVa[i])+" ")
        file.write("\n")
        for i in range(3):
            file.write(str(self.ACVaTFNumber)+" ")
        file.write("\n")  
            
        # Number of Eigenvalues computed
        if(self.AnalysisFlag ==3):
            file.write(str(self.Nev)+"\n")
        file.write("\n")  
        # Wing parameters
        # Number of keypoints
        file.write(str(nkp)+" ")
        # Number of members
        file.write(str(nmemb)+" ")
        # Number of condition points
        file.write("2 ")
        # Number of materials
        file.write(str(nmate)+" ")
        # Number of frame
        file.write(str(nframe)+" ")
        # Unused : number of condition members and number of distribution function
        file.write("0 0 ")
        # Number of timeFunction
        file.write(str(ntimefun)+" ")
        # Unused : number of curvatures
        file.write("0 \n\n")
        
        # Keypoint list
        for i in range(nkp):
            file.write(str(i+1)+" ")
            for j in range(3):
                file.write(str(self.Wing.GetKpList()[i][j])+" ")
            file.write("\n")     
        file.write("\n")                 
        
        # Member list
        for i in range(nmemb):
            # search for the mat and frame number
            imat = self.Wing.GetCrossSections().index(self.Wing.GetWingSections()[i].GetCrossSection())
            iframe = self.Wing.GetFrames().index(self.Wing.GetWingSections()[i].GetFrame())
            # Member number
            file.write(str(i+1)+" ")
            # first Kp number 
            file.write(str(i+1)+" ")
            # second Kp number 
            file.write(str(i+2)+" ")            
            
            
            # first material number 
            file.write(str(imat+1)+" ")
            # second material number 
            file.write(str(imat+1)+" ")
            # Frame number 
            file.write(str(iframe+1)+" ")
            # Number of divisions
            file.write(str(self.Wing.GetWingSections()[i].GetNumberOfElements())+" ")
            # Unused : Number of curvatures
            file.write("0\n")
        file.write("\n")    
        # Keypoints conditions
        file.write("1\n1 2 3 4 5 6\n0 0 0 0 0 0\n0 0 0 0 0 0\n0 0 0 0 0 0\n\n")
        file.write("{0}\n7 8 9 10 11 12\n0 0 0 0 0 0\n0 0 0 0 0 0\n0 0 0 0 0 0\n\n".format(nkp))
            
        # Cross section parameters
        for i in range(nmate):
            file.write(str(i+1)+"\n")
            # Flexibility Matrix
            FlexMat = self.Wing.GetCrossSections()[i].GetFlexibilityMatrix()
            for j in range(6):
                for k in range(6):
                    file.write(str(FlexMat[j][k])+" ")
                file.write("\n")
            file.write("\n")
            # Mass Matrix 
            MassMat = self.Wing.GetCrossSections()[i].GetMassMatrix()
            for j in range(6):
                for k in range(6):
                    file.write(str(MassMat[j][k])+" ")
                file.write("\n")
            # Aero parameters
            if (self.AeroFlag >0):
                # Vinf
                file.write(str(round(self.Vinf,8))+" ")
                # Rho
                file.write(str(self.Rho)+" ")
                # Chord
                file.write(str(self.Wing.GetWingSections()[i].GetChord())+" ")
                # Parameter a
                file.write(str(self.Wing.GetWingSections()[i].GetParameterA())+" ")
                # Alpha aircraft
                file.write(str(round(self.AlphaAC,8))+" ")
                # Beta aircraft
                file.write(str(round(self.BetaAC,8))+" ")
                # unused yet
                file.write("0\n")
                if self.GravFlag ==1:
                    file.write(str(round(self.Xcg,8))+"\n")
                file.write("\n")
                
        # Frame 
        for i in range (nframe):
            file.write(str(i+1)+"\n")
            Frame = self.Wing.GetFrames()[i].GetFrameMatrix()
            for j in range(3):
                for k in range(3):
                    file.write(str(Frame[j][k])+" ")
                file.write("\n")
            file.write("\n")    
        # Simu Time    
        if(self.AnalysisFlag <3):
            file.write(str(self.SimuStart)+" ")    
            file.write(str(self.SimuEnd)+"\n")    
            file.write("\n")            
        # TimeFunction
        if (ntimefun>0):
            for i in range(ntimefun):
                file.write(str(i+1)+"\n")
                TF = self.GetTimeFunction(i)
                nentrie = len(TF.GetFunctionEntries())
                file.write(str(TF.GetFunctionType())+"\n")
                file.write(str(TF.GetFunctionStart())+" ")
                file.write(str(TF.GetFunctionEnd())+"\n")
                file.write(str(nentrie)+"\n")                    
                for j in range(nentrie):
                    Entrie = TF.GetFunctionEntrie(j)
                    for k in range(len(Entrie)):
                        file.write(str(Entrie[k])+" ")
                            
                file.write("\n")                
        file.close()

    def RemoveInputFile(self):
        Command = "rm "+self.FileName+" input.ech "+self.FileName+".ini"
        sp.getoutput(Command)
        
    def WriteInitFile(self):
        # write the .ini filled with 0.0
        FileName = self.GetFileName()+".ini"
        file = open(FileName,"w")
        nmemb = len(self.Wing.GetWingSections())
        Ndiv = 0
        for i in range(nmemb):
            Ndiv = Ndiv + self.Wing.GetWingSections()[i].GetNumberOfElements()
        for i in range(2*Ndiv):
            file.write("0. 0. 0. 0. 0. 0.\n")
            
