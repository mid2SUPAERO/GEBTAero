# coding=UTF-8
import subprocess as sp

## This class is designed to write the input file .dat readable by the fortran computation code
class InputFile:
    def __init__(self,Name,AnalysisFlag,AeroFlag,GravFlag,Niter,Nstep,Nvtk,Nev,ACOmegaa,ACOmegaaTFNumber,ACVa,ACVaTFNumber,Wing,Vinf,Rho,AlphaAC,BetaAC,SimuStart,SimuEnd,Xcg=0.):
        ## the Name of the simulation
        self.Name = Name
        ## the name on the input file .dat
        self.FileName = Name+".dat"
        ## link to ioaero::analysis_flag
        self.AnalysisFlag = AnalysisFlag
        ## link to ioaero::aero_flag
        self.AeroFlag = AeroFlag
        ## link to ioaero::grav_flag
        self.GravFlag = GravFlag
        ## link to ioaero::niter
        self.Niter = Niter
        ## link to ioaero::nstep
        self.Nstep = Nstep
        ## link to ioaero::nvtk
        self.Nvtk = Nvtk
        ## link to ioaero::nev
        self.Nev = Nev
        ## link to ioaero::omega_a0
        self.ACOmegaa = ACOmegaa
        ## link to ioaero::omega_a_tf
        self.ACOmegaaTFNumber = ACOmegaaTFNumber
        ## link to ioaero::v_root_a0
        self.ACVa = ACVa
        ## link to ioaero::v_root_a_tf
        self.ACVaTFNumber = ACVaTFNumber
        ## the Wing used in the simulation
        self.Wing = Wing
        ## Upstream velocity of the simulation
        self.Vinf = Vinf
        ## air density of the simulation
        self.Rho = Rho
        ## Aircraft angle of attack
        self.AlphaAC = AlphaAC
        ## Aircraft yaw angle
        self.BetaAC = BetaAC
        ## Start time of the simulation
        self.SimuStart = SimuStart
        ## End time of the simulation
        self.SimuEnd = SimuEnd
        ## position of center of gravity in frame a
        self.Xcg = Xcg
        
        # Initialisation of time function list
        self.TimeFunctions = []
          
    def GetName(self):
        return self.Name
    
    def GetFileName(self):
        return self.FileName    
        
    def GetAnalysisFlag(self):
        return self.AnalysisFlag    
        
    ## Add a TimeFunction to the InputFile    
    def AppendTimeFunction(self,TimeFunction):
        self.TimeFunctions.append(TimeFunction)
    
    def GetTimeFunction(self,index):
        return self.TimeFunctions[index]
        
    ## Write a .dat input file readable by the fortran computation code    
    def WriteInputFile(self):
        nkp = len(self.Wing.GetKpList())
        nmemb = len(self.Wing.GetWingSections())
        nmate = len(self.Wing.GetCrossSections())
        nframe = len(self.Wing.GetFrames())
        ntimefun = len(self.TimeFunctions)
        file = open(self.FileName,"w")
        # first line of the solver
        file.write(str(int(self.AnalysisFlag)))
        file.write(" "+str(int(self.AeroFlag)))
        file.write(" "+str(int(self.GravFlag)))
        file.write(" "+str(int(self.Niter)))
        file.write(" "+str(int(self.Nstep)))
        file.write(" "+str(int(self.Nvtk))+"\n")
        # Omegaa and Va
        for i in range(3):
            file.write(str(float(self.ACOmegaa[i]))+" ")
        file.write("\n")
        for i in range(3):
            file.write(str(int(self.ACOmegaaTFNumber))+" ")
        file.write("\n")
        for i in range(3):
            file.write(str(float(self.ACVa[i]))+" ")
        file.write("\n")
        for i in range(3):
            file.write(str(int(self.ACVaTFNumber))+" ")
        file.write("\n")  
            
        # Number of Eigenvalues computed
        if(self.AnalysisFlag ==3):
            file.write(str(int(self.Nev))+"\n")
        file.write("\n")  
        # Wing parameters
        # Number of keypoints
        file.write(str(int(nkp))+" ")
        # Number of members
        file.write(str(int(nmemb))+" ")
        # Number of condition points
        file.write("2 ")
        # Number of materials
        file.write(str(int(nmate))+" ")
        # Number of frame
        file.write(str(int(nframe))+" ")
        # Unused : number of condition members and number of distribution function
        file.write("0 0 ")
        # Number of timeFunction
        file.write(str(int(ntimefun))+" ")
        # Unused : number of curvatures
        file.write("0 \n\n")
        
        # Keypoint list
        for i in range(nkp):
            file.write(str(i+1)+" ")
            for j in range(3):
                file.write(str(float(self.Wing.GetKpList()[i][j]))+" ")
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
            file.write(str(int(self.Wing.GetWingSections()[i].GetNumberOfElements()))+" ")
            # Unused : Number of curvatures
            file.write("0\n")
        file.write("\n")    
        # Keypoints conditions
        file.write("1\n1 2 3 4 5 6\n0 0 0 0 0 0\n0 0 0 0 0 0\n0 0 0 0 0 0\n\n")
        file.write("{0}\n7 8 9 10 11 12\n0 0 0 0 0 0\n0 0 0 0 0 0\n0 0 0 0 0 0\n\n".format(int(nkp)))
            
        # Cross section parameters
        for i in range(nmate):
            file.write(str(i+1)+"\n")
            # Flexibility Matrix
            FlexMat = self.Wing.GetCrossSections()[i].GetFlexibilityMatrix()
            for j in range(6):
                for k in range(6):
                    file.write(str(float(FlexMat[j][k]))+" ")
                file.write("\n")
            file.write("\n")
            # Mass Matrix 
            MassMat = self.Wing.GetCrossSections()[i].GetMassMatrix()
            for j in range(6):
                for k in range(6):
                    file.write(str(float(MassMat[j][k]))+" ")
                file.write("\n")
            # Aero parameters
            if (self.AeroFlag >0):
                # Vinf
                file.write(str(round(float(self.Vinf),8))+" ")
                # Rho
                file.write(str(float(self.Rho))+" ")
                # Chord
                file.write(str(float(self.Wing.GetWingSections()[i].GetChord()))+" ")
                # Parameter a
                file.write(str(float(self.Wing.GetWingSections()[i].GetParameterA()))+" ")
                # Alpha aircraft
                file.write(str(round(float(self.AlphaAC),8))+" ")
                # Beta aircraft
                file.write(str(round(float(self.BetaAC),8))+" ")
                # unused yet
                file.write("0\n")
                if self.GravFlag ==1:
                    file.write(str(round(float(self.Xcg),8))+"\n")
                file.write("\n")
                
        # Frame 
        for i in range (nframe):
            file.write(str(i+1)+"\n")
            Frame = self.Wing.GetFrames()[i].GetFrameMatrix()
            for j in range(3):
                for k in range(3):
                    file.write(str(float(Frame[j][k]))+" ")
                file.write("\n")
            file.write("\n")    
        # Simu Time    
        if(self.AnalysisFlag <3):
            file.write(str(float(self.SimuStart))+" ")    
            file.write(str(float(self.SimuEnd))+"\n")    
            file.write("\n")            
        # TimeFunction
        if (ntimefun>0):
            for i in range(ntimefun):
                file.write(str(i+1)+"\n")
                TF = self.GetTimeFunction(i)
                nentrie = len(TF.GetFunctionEntries())
                file.write(str(int(TF.GetFunctionType()))+"\n")
                file.write(str(float(TF.GetFunctionStart()))+" ")
                file.write(str(float(TF.GetFunctionEnd()))+"\n")
                file.write(str(nentrie)+"\n")                    
                for j in range(nentrie):
                    Entrie = TF.GetFunctionEntrie(j)
                    for k in range(len(Entrie)):
                        file.write(str(Entrie[k])+" ")
                            
                file.write("\n")                
        file.close()

    ## Remove the .dat input file, the .ini (initial condition for dynamic simulation) and the input.ech (summary of the input process)
    def RemoveInputFile(self):
        Command = "rm "+self.FileName+" input.ech "+self.FileName+".ini"
        sp.getoutput(Command)
        
    ## Write the .ini filled with 0.0 for a transient dynamic simulation with no initial displacement or velocity
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
            
