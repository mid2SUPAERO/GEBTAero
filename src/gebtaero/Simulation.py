# coding=UTF-8
import numpy as np
import subprocess as sp
import math
from .utils import *
from .InputFile import *
from .TimeFunction import *
from .GebtPlot import *
from numpy.linalg import norm as norm


class Simulation:
    """
    This class contains the methods design to find eigenvalues, static shape, flutter speed
    (modal or temporal), divergence speed of a Wing in a particular configuration (AoA, Vinf,...)
    """
    def __init__(self,Wing):
        self.Wing = Wing
        
    def GetWing(self):
        return self.Wing
        
    def Eigenvalues(self,Vinf,Rho,AlphaAC,BetaAC,AeroFlag,NumberOfModes,GravFlag=0,verbosity=0,arpack=0,vtk=0):
        Name = self.GetWing().GetName()+"Modes"
        VecNul = np.zeros([3])
        self.Input = InputFile(Name,3,AeroFlag,GravFlag,1000,100,vtk,NumberOfModes,VecNul,0,VecNul,0,self.Wing,Vinf,Rho,AlphaAC,BetaAC,0.,0.)       
        self.Input.WriteInputFile()
        FileName = self.Input.GetFileName()
        a = "arpack="+str(arpack)
        e = "eigenoutput=1"
        Command = ["gebtaero", "-p ",FileName,a,e]
        EigenValues= ReadModesFromPipe(Command,output="eigenvalues")
        EigenValues = np.array(EigenValues)
        self.Input.RemoveInputFile()    
        if verbosity == 1:
            print('[Modal damping / frequency]')
            print(EigenValues)
        return EigenValues

    def StaticLoads(self,Vinf,Rho,AlphaAC,BetaAC,GravFlag=0,verbosity=0):
        Name = self.GetWing().GetName()+"Static"+str(Vinf)+"ms"
        VecNul = np.zeros([3])
        self.Input = InputFile(Name,1,1,GravFlag,1000,1,0,0,VecNul,0,VecNul,0,self.Wing,Vinf,Rho,AlphaAC,BetaAC,0.,0.)       
        self.Input.WriteInputFile()
        FileName = self.Input.GetFileName()
        Command = ["gebtaero","-p",FileName]
        StaticLoads = ReadLoadsFromPipe(Command)
        Output = np.array(StaticLoads,dtype=float,copy=False).reshape(2,3)
        # ~ RemoveFiles()
        if verbosity == 1:
            print('The static loads in frame a is Fxa=',str(round(Output[0,0],4)),'N ; Fya=',str(round(Output[0,1],4)),'N ; Fza=',str(round(Output[0,2],4)),'N ; Mxa=',str(round(Output[1,0],4)),'N.m ; Mya=',str(round(Output[1,1],4)),'N.m ; Mza=',str(round(Output[1,2],4)),'N.m')
        self.Input.RemoveInputFile()
        return Output


    def ModalFlutterSpeed(self,Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,FreqLim,AlphaAC,BetaAC,KsiObj=1e-6,GravFlag=0,verbosity=0,arpack=1,ModesToCompute=20):
        Name = self.GetWing().GetName()+"ModalFlutter"
        VecNul = np.zeros([3])
        Velocity = Vmin
        NumberOfModes = ModesToCompute
        MinDamp = None
        self.Input = InputFile(Name,3,AeroFlag,GravFlag,1000,1,0,NumberOfModes,VecNul,0,VecNul,0,self.Wing,Velocity,Rho,AlphaAC,BetaAC,0.,0.)       
        self.Input.WriteInputFile()
        while (Vstep>DeltaV):
            if verbosity ==1:
                print("Velocity = {0}; Modes computed : {1}".format(Velocity,NumberOfModes) )
            v = "v="+str(Velocity)
            m = "m="+str(NumberOfModes)
            a = "arpack="+str(arpack)
            e = "eigenoutput=1"
            Command = ["gebtaero","-p",self.Input.GetFileName(),v,m,a,e]
            EigenValues = ReadModesFromPipe(Command,output="eigenvalues")
            EigenVal = np.array(EigenValues,dtype=float)    
            n = len(EigenVal)
            if n ==0:
                NumberOfModes = max(int(1.5*NumberOfModes),NumberOfModes+5)
                if verbosity ==1:
                    print("Modes to compute increased to {0}".format(NumberOfModes))
                continue # skip the end of the process to start again with the same velocity
            # adjusting the number of modes to calculate
            FreqMax = max(EigenVal[:,1])
            if FreqMax < FreqLim:
                NumberOfModes = max(int(1.4*NumberOfModes),NumberOfModes+5)
                if verbosity ==1:
                    print("Modes to compute increased to {0}".format(NumberOfModes))
                continue # skip the end of the process to start again with the same velocity
            elif FreqMax > 2* FreqLim:  
                NumberOfModes = int(0.9*NumberOfModes)   
            # reduced frequency calculation
            for i in range(n):
                if(math.isclose(EigenVal[i,1],0.)):
                    EigenVal[i,0] = 0.
                else:    
                    EigenVal[i,0] = EigenVal[i,0] /(-4*np.pi*abs(EigenVal[i,1]))      
      
            # search of the max damping
            FreqIndex = np.argwhere(EigenVal[:,1]>FreqLim)[0][0]    
            MinDamp = np.amin(EigenVal[:FreqIndex,0])
            # Determination of the new velocity depending on the damping
            if Velocity > Vmin:
                if (MinDamp>-KsiObj):

                    Velocity = Velocity+Vstep
                else:
                    Velocity = Velocity-0.5*Vstep
                    Vstep = 0.5*Vstep 
                if Velocity>Vmax :
                    if verbosity ==1:
                        print("no flutter detected in the velocity range")
                    return np.array([-1,-1])
            else:    
                Velocity = Velocity + Vstep        
            MinDampOld = MinDamp
            # ~ print(Velocity,MinDamp,Vstep,NumberOfModes)    
        self.Input.RemoveInputFile()
        # flutter frequency extraction
        FlutterIndex = np.argmin(EigenVal[:FreqIndex,0])
        FlutterFreq = EigenVal[FlutterIndex,1]   
        if verbosity ==1:
            print("The flutter speed is {0} ; the flutter frequency is {1} Hz converged with a tolerance of {2} m/s and a reduced damping of {3}".format(Velocity,FlutterFreq,DeltaV,MinDamp))
        return [Velocity,FlutterFreq]
        

    def ModalCriticalSpeed(self,Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=0,verbosity=0,mode=0):
        if mode ==3: #looking for both divergence and flutter speed
            CritVelocity1,CritFreq1 = self.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=GravFlag,verbosity=verbosity,mode=0)
            if CritFreq1 == 0.:# the first instability is a divergence
                CritVelocity2,CritFreq2 = self.ModalCriticalSpeed(Rho,CritVelocity1,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=GravFlag,verbosity=verbosity,mode=1)                
                return CritVelocity2,CritFreq2,CritVelocity1,CritFreq1 # the flutter is returned first
            else:# the first instability is a flutter
                CritVelocity2,CritFreq2 = self.ModalCriticalSpeed(Rho,CritVelocity1,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=GravFlag,verbosity=verbosity,mode=2)                
                return CritVelocity1,CritFreq1,CritVelocity2,CritFreq2 # the flutter is returned first

        Name = self.GetWing().GetName()+"ModalCritical"
        VecNul = np.zeros([3])
        Vstep = min(Vstep,0.5*(Vmax-Vmin))
        Velocity = Vmin+Vstep
        Vlist = []
        NumberOfModes = 1
        self.Input = InputFile(Name,3,AeroFlag,GravFlag,1000,1,0,NumberOfModes,VecNul,0,VecNul,0,self.Wing,Velocity,Rho,AlphaAC,BetaAC,0.,0.)       
        self.Input.WriteInputFile()
        while (Vstep>DeltaV and Velocity >=Vmin+DeltaV and Velocity <= Vmax):
            if verbosity ==1:
                print("Velocity = {0}; Modes computed = {1}".format(Velocity,NumberOfModes) )
            v = "v="+str(Velocity)
            m = "m="+str(NumberOfModes)
            a = "arpack=4"
            e = "eigenoutput=1"
            Command = ["gebtaero","-p",self.Input.GetFileName(),v,m,a,e]
            EigenValues = ReadModesFromPipe(Command,output="eigenvalues")
            EigenVal = np.array(EigenValues,dtype=float)
            # ~ if verbosity ==1:
                # ~ print(EigenVal)
            n = len(EigenVal)
            if (mode ==1): # looking for the first flutter mode
                # suppression of zero frequency mode
                i = 0
                while i<n:
                    if EigenVal[i,1] == 0.:
                        EigenVal = np.delete(EigenVal,i,axis=0)
                        n = n-1
                    else :
                        i = i+1
                n = len(EigenVal)
                if n ==0 :
                    NumberOfModes = max(int(1.5*NumberOfModes), NumberOfModes + 4)
                else:    
                    if (EigenVal[0,0]>0.):
                        CritDamp = EigenVal[0,0]   
                        CritFreq = EigenVal[0,1]
                        Vstep = 0.75*Vstep
                        Velocity = max(0.25*Vstep,Velocity - Vstep)
                    else:
                        if Velocity+Vstep in Vlist:
                            Vstep = 0.5*Vstep
                            Velocity = Velocity + Vstep
                        else:    
                            Velocity = Velocity + Vstep   
            elif (mode ==2): # looking for the first divergence mode
                if n ==0 :
                    NumberOfModes = max(int(1.5*NumberOfModes), NumberOfModes + 4)
                else:    
                    if (EigenVal[0,0]>0. and EigenVal[0,1]==0):  
                        CritDamp = EigenVal[0,0]   
                        CritFreq = EigenVal[0,1]
                        Vstep = 0.75*Vstep
                        Velocity = max(0.25*Vstep,Velocity - Vstep)
                    else:
                        if Velocity+Vstep in Vlist:
                            Vstep = 0.5*Vstep
                            Velocity = Velocity + Vstep
                        else:    
                            Velocity = Velocity + Vstep     
            else: # looking for the first unstable mode (divergence or flutter) mode 0
                if n ==0 :
                    NumberOfModes = max(int(1.5*NumberOfModes), NumberOfModes + 4)
                else:    
                    if (EigenVal[0,0]>0.):
                        CritDamp = EigenVal[0,0]      
                        CritFreq = EigenVal[0,1]
                        Vstep = 0.75*Vstep
                        Velocity = max(0.25*Vstep,Velocity - Vstep)
                    else:
                        if Velocity+Vstep in Vlist:
                            Vstep = 0.5*Vstep
                            Velocity = Velocity + Vstep
                        else:    
                            Velocity = Velocity + Vstep               
            # adaptation od the number of modes computed 
            if n >5:
                NumberOfModes = int(0.8*NumberOfModes)  
            Vlist.append(Velocity)
        
        self.Input.RemoveInputFile()
        if Velocity < Vmin+DeltaV:
            if verbosity ==1:
                print("the minimal speed is already unstable")
                Velocity = -1
            return[-1,-1]
        elif Velocity > Vmax:
            if verbosity ==1:
                print("the maximal speed is still stable")
                Velocity = -2    
            return[-2,-2]    
        elif mode ==2:
            if verbosity ==1:
                print("The divergence speed is {0} m/s computed with a precision of {1} m/s and a real part of {2}".format(round(Velocity,4),DeltaV,round(CritDamp,10)))
            return [Velocity,0.]
        elif mode == 1 : 
            if verbosity ==1:
                print("The Flutter speed is {0} m/s, the flutter frequency is {1} Hz computed with a precision of {2} m/s and a reduced damping of {3}".format(round(Velocity,4),round(CritFreq,4),DeltaV,round(CritDamp/(-4*np.pi*CritFreq),10)))
            return [Velocity,CritFreq]    
        else : 
            if verbosity ==1:
                print("The critical speed is {0} m/s, the critical frequency is {1} Hz computed with a precision of {2} m/s and a real part of {3}".format(round(Velocity,4),round(CritFreq,4),DeltaV,round(CritDamp,10)))
            return [Velocity,CritFreq]    
        


    def ModalDivergenceSpeed(self,Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,NormLim=1.,KsiObj=1e-6,GravFlag=0,verbosity=0):
        Name = self.GetWing().GetName()+"ModalFlutter"
        VecNul = np.zeros([3])
        Velocity = Vmin
        NumberOfModes = 20
        MinDamp = None
        self.Input = InputFile(Name,3,AeroFlag,GravFlag,1000,1,0,NumberOfModes,VecNul,0,VecNul,0,self.Wing,Velocity,Rho,AlphaAC,BetaAC,0.,0.)       
        self.Input.WriteInputFile()
        while (Vstep>DeltaV):
            if verbosity ==1:
                print("Velocity = {0}; Modes computed : {1}".format(Velocity,NumberOfModes) )
            v = "v="+str(Velocity)
            m = "m="+str(NumberOfModes)
            a = "arpack=2"
            e = "eigenoutput=1"
            Command = ["gebtaero","-p",self.Input.GetFileName(),v,m,a,e]
            EigenValues = ReadModesFromPipe(Command,output="eigenvalues")
            EigenVal = np.array(EigenValues,dtype=float)    
            n = len(EigenVal)
            if n ==0:
                NumberOfModes = max(int(1.2*NumberOfModes),NumberOfModes+5)
                if verbosity ==1:
                    print("Modes to compute increased to {0}".format(NumberOfModes))
                continue # skip the end of the process to start again with the same velocity
            # adjusting the number of modes to calculate
            NormMax = max(norm(EigenVal,axis=1))
            if NormMax < NormLim:
                NumberOfModes = max(int(1.2*NumberOfModes),NumberOfModes+5)
                if verbosity ==1:
                    print("Modes to compute increased to {0}".format(NumberOfModes))
                continue # skip the end of the process to start again with the same velocity
            elif NormMax > 2* NormLim:  
                NumberOfModes = int(0.9*NumberOfModes)         
            # search of the minimal damping
            Ind = np.argwhere(EigenVal[:,1] ==0.).ravel()
            if len(Ind) ==0:
                MinDamp = -1
            else :
                MinDamp = np.amax(EigenVal[Ind,0])
            # Determination of the new velocity depending on the damping
            if Velocity > Vmin:
                if MinDamp<KsiObj :
                    Velocity = Velocity+Vstep
                else:
                    Velocity = Velocity-0.5*Vstep
                    Vstep = 0.5*Vstep 
                if Velocity>Vmax :
                    if verbosity ==1:
                        print("no divergence detected in the velocity range")
                    return np.array([-1,-1])
            else:    
                Velocity = Velocity + Vstep        
        self.Input.RemoveInputFile()
        if verbosity ==1:
            print("The divergence speed is {0} converged with a tolerance of {1} m/s and a real part of {2}".format(Velocity,DeltaV,MinDamp))
        return Velocity
        
        
        
    def ModalFlutterSpeedSorted(self,Rho,Vmax,DeltaV,AeroFlag,ModesToCompute,ModesToPlot,AlphaAC,BetaAC,KsiObj=1e-6,CorrCoef=0.9,GravFlag=0,verbosity=0,arpack=0):
        if AeroFlag == 3:
            ModesToCompute = max(ModesToCompute,5*ModesToPlot)
        else:    
            ModesToCompute = max(ModesToCompute,3*ModesToPlot)
        Name = self.GetWing().GetName()+"ModalFlutter"
        VecNul = np.zeros([3])
        self.Input = InputFile(Name,3,AeroFlag,GravFlag,1000,1,0,ModesToCompute,VecNul,0,VecNul,0,self.Wing,1.,Rho,AlphaAC,BetaAC,0.,0.)       
        self.Input.WriteInputFile()
        VstepMax = min(Vmax/10.,10.)
        Vstep = VstepMax
        VstepCounter = 0
        Vinf = 0.
        corr_flag = False
        Velocity = []
        Modes = []
        while Vstep>DeltaV and Vinf<Vmax:
            v = "v="+str(Vinf)
            m = "m="+str(ModesToCompute)
            a = "arpack="+str(arpack)
            Command = ["gebtaero","-p",self.Input.GetFileName(),v,m,a]
            EigenModes,errcode = ReadModesFromPipe(Command,output="modes")
            if errcode ==106:
                ModesToCompute = max(int(1.2*ModesToCompute),ModesToCompute+5)
                if verbosity ==1:
                    print("ModesToCompute increased to {0}; gebtaero err 106".format(ModesToCompute))
                if ModesToCompute > 100*ModesToPlot:
                    raise RuntimeError("Unable to perform the eigenvalue calculation, the number of beam element is probably too small")
                else:    
                    return self.ModalFlutterSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToCompute,ModesToPlot,AlphaAC,BetaAC,KsiObj=KsiObj,CorrCoef=CorrCoef,GravFlag=GravFlag,verbosity=verbosity,arpack=arpack)
            n = len(EigenModes)
            if n == 0:
                ModesToCompute = max(int(1.2*ModesToCompute),ModesToCompute+5)
                continue    
            EigenValues= []
            for i in range(n):
                EigenValues.append(EigenModes[i][0])
            EigenVal = np.array(EigenValues,dtype=float)
            #initialisation for Vinf = Vmin:
            if Vinf == 0. :
                # selecting the first ModesToPlot with lowest and nonzero frequency
                IndexSort = np.argsort(abs(EigenVal[:,1]))
                MinInd = int(np.argwhere(EigenVal[IndexSort,1]>0)[0])
                IndexOld = IndexSort[MinInd:MinInd+ModesToPlot]
                # test if ModesToCompute is sufficient
                if len(IndexOld)<ModesToPlot:
                    ModesToCompute = max(int(1.2*ModesToCompute),ModesToCompute+5)
                    if verbosity ==1:
                        print("ModesToCompute increased to {0} (n too small)".format(ModesToCompute))
                    return self.ModalFlutterSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToCompute,ModesToPlot,AlphaAC,BetaAC,KsiObj=KsiObj,CorrCoef=CorrCoef,GravFlag=GravFlag,verbosity=verbosity,arpack=arpack)
                EigenModesOld = EigenModes
                EigenValOld = EigenVal
                VinfOld = Vinf
                Vinf = Vinf+Vstep
                Velocity.append(VinfOld)
                Modes.append(EigenVal[IndexOld,:])
            elif(max(norm(EigenVal,axis=1))<max(norm(EigenValOld[IndexOld,:],axis=1))):
                ModesToCompute = max(int(1.2*ModesToCompute),ModesToCompute+5)
            else:    
                IndexNew = []
                Corr = CorrelateTab(EigenModesOld,EigenModes,IndexOld)
                corr_flag = False
                for i in range(ModesToPlot): 
                    CorrLine = Corr[IndexOld[i],:,0]*Corr[IndexOld[i],:,1]                                              
                    CorrMax = CorrLine.argmax()
                    CorrVec = Corr[IndexOld[i],CorrMax,0]
                    CorrVal = Corr[IndexOld[i],CorrMax,1]
                    if CorrLine[CorrMax]>CorrCoef:
                        IndexNew.append(CorrMax)
                    elif (max(norm(EigenVal,axis=1))>2*max(norm(EigenValOld[IndexOld,:],axis=1))) :
                        Vstep = 0.5*Vstep
                        Vinf = Vinf - Vstep
                        corr_flag = True
                        if verbosity ==1:
                            print("Vstep decreased to {0}".format(Vstep))
                        break    
                    else:
                        ModesToCompute = max(int(1.2*ModesToCompute),ModesToCompute+5)
                        corr_flag = True
                        if verbosity ==1:
                            # ~ print(max(norm(EigenVal,axis=1))<2*max(norm(EigenValOld[IndexOld,:],axis=1)))
                            print("ModesToCompute increased to {0}".format(ModesToCompute))
                        break    
                if corr_flag:
                    continue        
                # reduced frequency calculation
                for i in range(ModesToPlot):
                    if(math.isclose(EigenVal[IndexNew[i],1],0.)):
                        EigenVal[IndexNew[i],0] = 0.
                    else:    
                        EigenVal[IndexNew[i],0] = EigenVal[IndexNew[i],0] /(-4*np.pi*abs(EigenVal[IndexNew[i],1]))       
                Velocity.append(Vinf)
                Modes.append(EigenVal[IndexNew,:])
                #minimal damping
                MinDamp = min(EigenVal[IndexNew,0])
                if MinDamp <-KsiObj: #Vinf unstable
                    Vstep = 0.5*Vstep
                    Vinf = Vinf - Vstep  
                else: #Vinf stable
                    Vinf = Vinf + Vstep
                    VstepCounter = VstepCounter + 1
                    if VstepCounter >4:
                        Vstep = min(VstepMax,2*Vstep)
                        VstepCounter = 0
                    if max(norm(EigenVal,axis=1))>3*max(norm(EigenValOld[IndexOld,:],axis=1)):
                        ModesToCompute = int(0.9*ModesToCompute)
                #update old values
                EigenModesOld = EigenModes
                EigenValOld = EigenVal
                VinfOld = Vinf 
                IndexOld=IndexNew   
            if verbosity ==1:
                print(VinfOld,ModesToCompute)  
        if corr_flag:
            raise RuntimeError("unable to follow the modes {0} from Vinf = {1} to {2} m/s".format(EigenValOld[IndexOld[i],:],Velocity[-1],Vinf))    
        if Vinf>=Vmax:
            return -2,-2
        Freq =  abs(EigenVal[IndexNew[EigenVal[IndexNew,0].argmin()],1])                           
        if verbosity == 1:
            print('The flutter speed is ',str(round(Vinf,2)),' ; the flutter frequency is ',str(round(Freq,2)),' Hz converged with a tolerance of ',str(round(DeltaV,4)),'m/s and a reduced damping of ',str(MinDamp))
        for i in range (len(Modes)):
            for j in range(len(Modes[i])):
                if Modes[i][j][0] == 0.:
                    Modes[i][j][0] = None
        Modes = np.array(Modes)
        return Vinf,Freq,[Velocity,Modes] 
   
   
    def ModalDivergenceSpeedSorted(self,Rho,Vmax,DeltaV,AeroFlag,ModesToCompute,ModesToPlot,AlphaAC,BetaAC,CorrCoef=0.95,GravFlag=0,verbosity=0):
        if AeroFlag == 3:
            ModesToCompute = max(ModesToCompute,5*ModesToPlot)
        else:    
            ModesToCompute = max(ModesToCompute,3*ModesToPlot)
        Name = self.GetWing().GetName()+"ModalFlutter"
        VecNul = np.zeros([3])
        self.Input = InputFile(Name,3,AeroFlag,GravFlag,1000,1,0,ModesToCompute,VecNul,0,VecNul,0,self.Wing,1.,Rho,AlphaAC,BetaAC,0.,0.)       
        self.Input.WriteInputFile()
        VstepMax = min(Vmax/10.,10.)
        Vstep = VstepMax
        VstepCounter = 0
        Vinf = 0.
        corr_flag = False
        Velocity = []
        Modes = []
        while Vstep>DeltaV and Vinf<Vmax:
            v = "v="+str(Vinf)
            m = "m="+str(ModesToCompute)
            Command = ["gebtaero","-p",self.Input.GetFileName(),v,m]
            EigenModes,errcode = ReadModesFromPipe(Command,output="modes")
            if errcode ==106:
                ModesToCompute = max(int(1.2*ModesToCompute),ModesToCompute+5)
                print("ModesToCompute increased to {0}".format(ModesToCompute+2))
                return self.ModalDivergenceSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToCompute,ModesToPlot,AlphaAC,BetaAC,CorrCoef=CorrCoef,GravFlag=GravFlag,verbosity=verbosity)
            n = len(EigenModes)
            EigenValues= []
            for i in range(n):
                EigenValues.append(EigenModes[i][0])
            EigenVal = np.array(EigenValues,dtype=float)
            #initialisation for Vinf = Vmin:
            if Vinf == 0. :
                # selecting the first ModesToPlot with lowest and nonzero frequency
                IndexSort = np.argsort(abs(EigenVal[:,1]))
                MinInd = int(np.argwhere(EigenVal[IndexSort,1]>0)[0])
                IndexOld = IndexSort[MinInd:MinInd+ModesToPlot]  
                # test if ModesToCompute is sufficient
                if len(IndexOld)<ModesToPlot:
                    ModesToCompute = max(int(1.2*ModesToCompute),ModesToCompute+5)
                    print("ModesToCompute increased to {0}".format(ModesToCompute))
                    return self.ModalDivergenceSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToCompute,ModesToPlot,AlphaAC,BetaAC,CorrCoef=CorrCoef,GravFlag=GravFlag,verbosity=verbosity)
                EigenModesOld = EigenModes
                EigenValOld = EigenVal
                VinfOld = Vinf
                Vinf = Vinf+Vstep
                Velocity.append(VinfOld)
                Modes.append(EigenVal[IndexOld,:])
            elif(max(norm(EigenVal,axis=1))<max(norm(EigenValOld[IndexOld,:],axis=1))):
                ModesToCompute = max(int(1.2*ModesToCompute),ModesToCompute+5)
            else:    
                IndexNew = []
                Corr = CorrelateTab(EigenModesOld,EigenModes,IndexOld)
                corr_flag = False
                for i in range(ModesToPlot): 
                    CorrLine = Corr[IndexOld[i],:,0]*Corr[IndexOld[i],:,1]                        
                    # ~ CorrLine = Corr[IndexOld[i],:,1]                        
                    CorrMax = CorrLine.argmax()
                    CorrVec = Corr[IndexOld[i],CorrMax,0]
                    CorrVal = Corr[IndexOld[i],CorrMax,1]
                    if CorrLine[CorrMax]>CorrCoef:
                        IndexNew.append(CorrMax)
                    elif (max(norm(EigenVal,axis=1))>2*max(norm(EigenValOld[IndexOld,:],axis=1))) :
                        Vstep = 0.5*Vstep
                        Vinf = Vinf - Vstep
                        corr_flag = True
                        if verbosity ==1:
                            print("Vstep decreased to {0}".format(Vstep))
                        break    
                    else:
                        ModesToCompute = max(int(1.2*ModesToCompute),ModesToCompute+5)
                        corr_flag = True
                        if verbosity ==1:
                            print("ModesToCompute increased to {0}".format(ModesToCompute))
                        break    
                if corr_flag:
                    continue        
                Velocity.append(Vinf)
                Modes.append(EigenVal[IndexNew,:])
                #maximal Real part
                MaxRe = -1
                for i in range(ModesToPlot):
                    if (math.isclose(EigenVal[IndexNew[i],1],0.)):
                        MaxRe = max(MaxRe,EigenVal[IndexNew[i],0])
                if MaxRe >0.: #Vinf unstable
                    Vstep = 0.5*Vstep
                    Vinf = Vinf - Vstep  
                else: #Vinf stable
                    Vinf = Vinf + Vstep
                    VstepCounter = VstepCounter + 1
                    if VstepCounter >4:
                        Vstep = min(VstepMax,2*Vstep)
                        VstepCounter = 0
                    if max(norm(EigenVal,axis=1))>2*max(norm(EigenValOld[IndexOld,:],axis=1)):
                        ModesToCompute = int(0.9*ModesToCompute)
                #update old values
                EigenModesOld = EigenModes
                EigenValOld = EigenVal
                VinfOld = Vinf 
                IndexOld=IndexNew   
            if verbosity ==1:
                print(VinfOld)  
        if corr_flag:
            raise RuntimeError("unable to follow the modes at Vinf = {0}".format(Vinf))    
        if Vinf>=Vmax:
            return -2,-2                     
        if verbosity == 1:
            print('The divergence speed is ',str(round(Vinf,2)),' converged with a tolerance of ',str(round(DeltaV,4)))
        for i in range (len(Modes)):
            for j in range(len(Modes[i])):
                if Modes[i][j][0] == 0.:
                    Modes[i][j][0] = None
        Modes = np.array(Modes)
        self.Input.RemoveInputFile()
        return Vinf,[Velocity,Modes] 
    
    def EigenTabSorted(self,Rho,Vmax,DeltaV,Nstep,AeroFlag,ModesToCompute,ModesToPlot,AlphaAC,BetaAC,CorrCoef=0.9,GravFlag=0,verbosity=0):
        if AeroFlag == 3:
            ModesToCompute = max(ModesToCompute,5*ModesToPlot)
        else:    
            ModesToCompute = max(ModesToCompute,3*ModesToPlot)
        Name = self.GetWing().GetName()+"ModalPlot"
        VecNul = np.zeros([3])
        self.Input = InputFile(Name,3,AeroFlag,GravFlag,1000,1,0,ModesToCompute,VecNul,0,VecNul,0,self.Wing,1.,Rho,AlphaAC,BetaAC,0.,0.)       
        self.Input.WriteInputFile()
        VstepMax = float(Vmax/Nstep)
        Vstep = VstepMax
        VstepCounter = 0
        Vinf = 0.
        corr_flag = False
        Velocity = []
        Modes = []
        while Vinf<Vmax+0.5*Vstep and Vstep>DeltaV:
            v = "v="+str(Vinf)
            m = "m="+str(ModesToCompute)
            a = "arpack="+str(2)
            Command = ["gebtaero","-p",self.Input.GetFileName(),v,m,a]
            EigenModes,errcode = ReadModesFromPipe(Command,output="modes")
            if errcode ==106:
                ModesToCompute = max(int(1.2*ModesToCompute),ModesToCompute+5)
                if verbosity ==1:
                    print("ModesToCompute increased to {0}; gebtaero err 106".format(ModesToCompute))
                if ModesToCompute > 100*ModesToPlot:
                    raise RuntimeError("Unable to perform the eigenvalue calculation, the number of beam element is probably too small")
                else:    
                    return self.EigenTabSorted(Rho,Vmax,Nstep,AeroFlag,ModesToCompute,ModesToPlot,AlphaAC,BetaAC,CorrCoef=CorrCoef,GravFlag=GravFlag,verbosity=verbosity)
            n = len(EigenModes)
            if n < ModesToPlot:
                ModesToCompute = max(int(1.2*ModesToCompute),ModesToCompute+5)
                continue    
            EigenValues= []
            ModeToAdd = -1
            ModeToRemove = -1
            for i in range(n):
                EigenValues.append(EigenModes[i][0])
            EigenVal = np.array(EigenValues,dtype=float)
            #initialisation for Vinf = Vmin:
            if Vinf == 0. :
                # determining the first ModesToPlot aerolastic modes to follow (the non zero lowest frequency modes)
                IndexSort = np.argsort(abs(EigenVal[:,1]))
                MinInd = int(np.argwhere(EigenVal[IndexSort,1]>0)[0])
                IndexOld = IndexSort[MinInd:MinInd+ModesToPlot]
                # backup the values
                EigenModesOld = EigenModes
                EigenValOld = EigenVal
                VinfOld = Vinf
                Vinf = Vinf+Vstep
                Velocity.append(VinfOld)
                Modes.append(EigenVal[IndexOld,:])
            elif(max(norm(EigenVal,axis=1))<max(norm(EigenValOld[IndexOld,:],axis=1))): # not enough modes computed
                ModesToCompute = max(int(1.2*ModesToCompute),ModesToCompute+5)
            else: # general case
                IndexNew = [] # array to store the index of the followed modes of the next Velocity step
                Corr = CorrelateTab(EigenModesOld,EigenModes,IndexOld)
                corr_flag = False
                for i in range(ModesToPlot): 
                    CorrLine = Corr[IndexOld[i],:,0]*Corr[IndexOld[i],:,1]                  
                    # ~ CorrLine = Corr[IndexOld[i],:,1]                        
                    CorrMax = CorrLine.argmax()
                    CorrVec = Corr[IndexOld[i],CorrMax,0]
                    CorrVal = Corr[IndexOld[i],CorrMax,1]
                    if CorrLine[CorrMax]>CorrCoef :
                        IndexNew.append(CorrMax)
                        # ajout du suivi d'un deuxième mode réel dans le cas d'un mode complexe de fréquence qui s'annule
                        if EigenValOld[IndexOld[i],1] > 0. and EigenVal[CorrMax,1] == 0.:
                            CorrLine[CorrMax] = 0.
                            ModeToAdd = CorrLine.argmax()   
                            Corr[:,CorrMax,:]=0.
                        # fusion de 2 modes réels en 1 complexe
                        # ~ elif EigenValOld[IndexOld[i],1] == 0. and EigenVal[CorrMax,1] > 0.:   
                            # ~ print(i,IndexNew[i],EigenValOld[IndexOld[i],:],EigenVal[CorrMax,:])
                            # ~ print(IndexNew)
                        # ~ else :
                            # ~ Corr[:,CorrMax,:]=0.
                    elif Vstep<2*DeltaV :    
                        IndexNew.append(CorrMax)
                    elif (max(norm(EigenVal,axis=1))>2*max(norm(EigenValOld[IndexOld,:],axis=1))):
                        Vstep = 0.5*Vstep
                        Vinf = Vinf - Vstep
                        corr_flag = True
                        if verbosity ==1:
                            print("Vstep decreased to {0}".format(Vstep))
                        break    
                    else:
                        ModesToCompute = max(int(1.2*ModesToCompute),ModesToCompute+5)
                        corr_flag = True
                        if verbosity ==1:
                            print("ModesToCompute increased to {0}".format(ModesToCompute))
                        break  
                # adding the supplementary real mode
                if ModeToAdd > -1:
                    ModesToPlot = ModesToPlot +1   
                    IndexNew.append(ModeToAdd)   
                if corr_flag:
                    continue        
                # ~ # reduced frequency calculation
                # ~ for i in range(ModesToPlot):
                    # ~ if(math.isclose(EigenVal[IndexNew[i],1],0.)):
                        # ~ EigenVal[IndexNew[i],0] = None
                    # ~ else:    
                        # ~ EigenVal[IndexNew[i],0] = EigenVal[IndexNew[i],0] /(-4*np.pi*abs(EigenVal[IndexNew[i],1]))       
                Velocity.append(Vinf)
                Modes.append(EigenVal[IndexNew,:])
                if VstepCounter >4:
                    Vstep = min(VstepMax,2*Vstep)
                    VstepCounter = 0
                Vinf = Vinf+Vstep    
                VstepCounter = VstepCounter+1
                if max(norm(EigenVal,axis=1))>2*max(norm(EigenValOld[IndexOld,:],axis=1)):
                        ModesToCompute = int(0.9*ModesToCompute)
                #update old values
                EigenModesOld = EigenModes
                EigenValOld = EigenVal
                VinfOld = Vinf 
                IndexOld=IndexNew   
                # ~ print(IndexNew)
            if verbosity ==1:
                print("Velocity = {0}; {1} modes computed".format(Vinf,ModesToCompute))  
        if corr_flag:
            raise RuntimeError("unable to follow the modes at Vinf = {0}".format(Vinf))    
        # completion of Modes before transform it into np array
        for i in range(len(Velocity)):
            for j in range(ModesToPlot-len(Modes[i])):
                Modes[i]=np.append(Modes[i],[[None,None]],axis=0)
        Modes = np.array(Modes) 
        return [Velocity,Modes]
            
    def EigenTab(self,Rho,Vmin,Vmax,Nstep,AeroFlag,ModesToPlot,ModesToCompute,AlphaAC,BetaAC,GravFlag=0,verbosity=0):
        if AeroFlag == 3:
            ModesToCompute = max(ModesToCompute,2*ModesToPlot)
        else:    
            ModesToCompute = max(ModesToCompute,2*ModesToPlot) 
        Name = self.GetWing().GetName()+"ModalPlot"
        VecNul = np.zeros([3])
        self.Input = InputFile(Name,3,AeroFlag,GravFlag,1000,1,0,ModesToCompute,VecNul,0,VecNul,0,self.Wing,1.,Rho,AlphaAC,BetaAC,0.,0.)       
        self.Input.WriteInputFile()
        VstepMax = float(Vmax/Nstep)
        Vstep = VstepMax
        Vinf = Vmin
        Velocity = []
        Modes = []
        while Vinf<=Vmax :
            v = "v="+str(Vinf)
            m = "m="+str(ModesToCompute)
            Command = ["gebtaero","-p",self.Input.GetFileName(),v,m]
            EigenModes = ReadModesFromPipe(Command,output="eigenvalues")
            # ~ if errcode ==106:
                # ~ ModesToCompute = max(int(1.2*ModesToCompute),ModesToCompute+5)
                # ~ if verbosity ==1:
                    # ~ print("ModesToCompute increased to {0}; gebtaero err 106".format(ModesToCompute))
                # ~ if ModesToCompute > 100*ModesToPlot:
                    # ~ raise RuntimeError("Unable to perform the eigenvalue calculation, the number of beam element is probably too small")
                # ~ else:    
                    # ~ return self.EigenTab(Rho,Vmin,Vmax,Nstep,AeroFlag,ModesToPlot,ModesToCompute,AlphaAC,BetaAC,GravFlag=GravFlag,verbosity=verbosity)
            n = len(EigenModes)
            if n < ModesToPlot:
                ModesToCompute = max(int(1.2*ModesToCompute),ModesToCompute+5)
                continue    
            EigenValues= []
            EigenVal = np.array(EigenModes,dtype=float)
            Velocity.append(Vinf)
            Modes.append(EigenVal)
            Vinf = Vinf+Vstep    
            if n>3*ModesToPlot:
                ModesToCompute = int(0.9*ModesToCompute)
            if verbosity ==1:
                print("Velocity = {0}; {1} modes computed".format(Vinf,ModesToCompute))  
        # ~ Modes = np.array(Modes)
        return [Velocity,Modes]
            
            
    def DeformedModalFlutterSpeed(self,Rho,Vmin,Vmax,DeltaV,AeroFlag,FreqLim,NumberOfModes,Ksitol,Lifttol,BetaAC,verbosity=0):
        Weight = self.Wing.GetWeight()
        Surface = self.Wing.GetSurface()
        Vstep = 0.1*(Vmax-Vmin)
        # Initialisation with the undeformed flutter speed
        VFlutter = self.ModalFlutterSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,FreqLim,0.,BetaAC,Ksitol)[0]
        VelocityOld = VFlutter
        # Equilibrium AlphaAC 
        AlphaEq = self.EquilibriumAoA(Rho,VFlutter,BetaAC,Lifttol)
        VFlutter = self.ModalFlutterSpeed(Rho,Vmin,Vmax,DeltaV,AeroFlag,NumberOfModes,AlphaEq,BetaAC,Ksitol)[0]
        counter = 0
        while(abs(VFlutter-VelocityOld)>DeltaV):
            counter = counter+1
            if (counter >100):
                raise RuntimeError('the algorithm is unable to converge : more than ',counter,' iterations yet')
            VelocityOld = VFlutter
            AlphaEq = self.EquilibriumAoA(Rho,VFlutter,BetaAC,Lifttol)
            Flutter = self.ModalFlutterSpeed(Rho,Vmin,Vmax,DeltaV,AeroFlag,NumberOfModes,AlphaEq,BetaAC,Ksitol)
            VFlutter = Flutter[0]
        if verbosity==1 :
            print('Flutter speed = ',str(round(VFlutter,2)),' m/s ; Flutter frequency = ',str(round(Flutter[1],2))+' Hz ; Equilibirum AoA = ',str(round(AlphaEq,4)),' deg')
        return Flutter
            
    def DeformedModalFlutterSpeedSorted(self,Rho,Vmax,DeltaV,AeroFlag,ModesToCompute,ModesToPlot,Lifttol,BetaAC,KsiObj=1e-6,verbosity=0):
        Weight = self.Wing.GetWeight()
        Surface = self.Wing.GetSurface()
        # Initialisation with the undeformed flutter speed
        VFlutter = self.ModalFlutterSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToCompute,ModesToPlot,0.,BetaAC,KsiObj)[0]
        VelocityOld = VFlutter
        # Equilibrium AlphaAC 
        AlphaEq = self.EquilibriumAoA(Rho,VFlutter,BetaAC,Lifttol)
        VFlutter = self.ModalFlutterSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToCompute,ModesToPlot,AlphaEq,BetaAC)[0]
        counter = 0
        while(abs(VFlutter-VelocityOld)>DeltaV):
            counter = counter+1
            if (counter >100):
                raise RuntimeError('the algorithm is unable to converge : more than ',counter,' iterations yet')
            VelocityOld = VFlutter
            AlphaEq = self.EquilibriumAoA(Rho,VFlutter,BetaAC,Lifttol)
            Flutter = self.ModalFlutterSpeedSorted(Rho,Vmax,DeltaV,AeroFlag,ModesToCompute,ModesToPlot,AlphaEq,BetaAC)
            VFlutter = Flutter[0]
        if verbosity==1 :
            print('Flutter speed = ',str(round(VFlutter,2)),' m/s ; Flutter frequency = ',str(round(Flutter[1],2))+' Hz ; Equilibirum AoA = ',str(round(AlphaEq,4)),' deg')
        return Flutter
        
        
    def EquilibriumAoA(self,Rho,Vinf,BetaAC,tolerance,verbosity=0):
        Weight = self.GetWing().GetWeight()
        Surface = self.GetWing().GetSurface()
        # Initialisation
        if (math.isclose(Vinf,0.)):
            return 0.
        else:
            AlphaEq = Weight/(Rho*np.pi*Surface*Vinf**2)
            Lift = -self.StaticLoads(Vinf,Rho,AlphaEq,BetaAC)[0][2]
            
            while((Lift-Weight)/max(Weight,1e-3)>tolerance):
                AlphaEq = Weight/max(Lift,1e-3)*AlphaEq
                Lift = -self.StaticLoads(Vinf,Rho,AlphaEq,BetaAC)[0][2]
        if verbosity == 1:
            print('Equilibirum Angle of Attack = ',str(round(AlphaEq,4)),' deg for a flow velocity of ',str(round(Vinf,2)),' corresponding to a Lift of ',str(round(Lift,2)))        
        return AlphaEq
  
        
    def TemporalFlutterSpeed(self,Rho,Vmin,Vmax,DeltaV,AeroFlag,ModesToPlot,AlphaAC,BetaAC,Ksitol,NbPeriod,StepByPeriod,CoefPerturb,GravFlag=0,verbosity=0):  
        # Determination of the flutter speed using the transient dynamic capacity of gebtaero
        Vstep = 0.1*(Vmax-Vmin)
        ModalFlutter = self.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=GravFlag,verbosity=verbosity,mode=1)
        VFlutter = round(ModalFlutter[0],4)
        Vinf = round(1.1*VFlutter,4)
        Vstep = round(0.1*VFlutter,4)
        FreqFlutter = round(ModalFlutter[1],4)
        if (math.isclose(FreqFlutter,0.)):
            raise RuntimeError('the modal flutter frequency is equal to zero')
        else:
            PeriodFlutter = round(1./FreqFlutter,4)
        #simulation time
        Time = NbPeriod*PeriodFlutter
        DeltaT = float(PeriodFlutter)/StepByPeriod        
        #initial perturbation parameters
        VecNul = np.zeros([3])
        Vz = round(CoefPerturb*Vinf,8)
        #~ Vz = CoefPerturb*Vinf
        #time simulation
        Nstep = int(float(Time)/DeltaT)
        flutter_flag = 1
        # Time Function.
        TF = TimeFunction(1,0.,PeriodFlutter)
        TF.AppendFunctionEntrieHarmonic(1.,PeriodFlutter,0.)
        # Input file creation
        Name = self.GetWing().GetName()+"TempFlutter"
        self.Input = InputFile(Name,2,AeroFlag,0,1000,Nstep,0,0,VecNul,1,VecNul,1,self.Wing,Vinf,Rho,AlphaAC,BetaAC,0.,1.)       
        self.Input.AppendTimeFunction(TF)
        self.Input.WriteInputFile()
        self.Input.WriteInitFile()
        
        while (flutter_flag == 1 and Vinf>0) :
            # call of gebtaero, the velocity of the input file is override by v= of the command
            if verbosity ==1:
                print('temporal simulation at Vinf=',str(Vinf),' ; number of time steps=',str(Nstep))
            Command = "gebtaero -s "+self.Input.GetFileName()+" v="+str(Vinf)+" aero="+str(AeroFlag)+" time="+str(Time)+" vz="+str(Vz)+" nstep="+str(Nstep)+" tfe="+str(PeriodFlutter)+" tfper="+str(PeriodFlutter)
            Output = sp.getoutput(Command)
            Output = Output.split()
            if (Output[0] == "no"):
                if (Vstep > DeltaV):
                    Vinf = Vinf + 0.75*Vstep
                    Vstep = max(DeltaV,0.5*Vstep)
                    Vinf = Vinf - Vstep
                    Vz =CoefPerturb*Vinf
                else:
                    flutter_flag = 0
                    break
            else:
                Vinf = Vinf - 0.25*Vstep
                Vz = CoefPerturb*Vinf
            Vinf = round(Vinf,4)    
            Vz = round(Vz,6)
        self.Input.RemoveInputFile()     
        return(Vinf)

    def FlutterVtk(self,Rho,Vmin,Vmax,DeltaV,AeroFlag,AlphaAC,BetaAC,NbPeriod,StepByPeriod,CoefPerturb,CoefVinf,Nvtk,FlutterLimit=0.4,GravFlag=0,verbosity=0):  
        # Computation of vtk files corresponding to a temporal flutter
        import shutil
        Vstep = 0.1*(Vmax-Vmin)
        ModalFlutter = self.ModalCriticalSpeed(Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,AlphaAC,BetaAC,GravFlag=GravFlag,verbosity=verbosity,mode=1)
        VFlutter = round(ModalFlutter[0],4)
        Vinf = round(CoefVinf*VFlutter,4)
        FreqFlutter = round(ModalFlutter[1],4)
        if (math.isclose(FreqFlutter,0.)):
            raise RuntimeError('the modal flutter frequency is equal to zero')
        else:
            PeriodFlutter = round(1./FreqFlutter,4)
        #simulation time
        Time = NbPeriod*PeriodFlutter
        DeltaT = float(PeriodFlutter)/StepByPeriod        
        #initial perturbation parameters
        VecNul = np.zeros([3])
        #~ Vz = round(CoefPerturb*Vinf,4)
        Vz = CoefPerturb*Vinf
        #time simulation
        Nstep = int(float(Time)/DeltaT)
        flutter_flag = 1
        # Time Function.
        TF = TimeFunction(1,0.,PeriodFlutter)
        TF.AppendFunctionEntrieHarmonic(1.,PeriodFlutter,0.)
        # Name creation
        Name = self.GetWing().GetName()+"TempFlutter"
        # suppression of previous vtk folder
        if (Nvtk > 0):
            vtk_folder = Name+".dat_vtk"
            shutil.rmtree(vtk_folder,True)
        # Input file creation
        self.Input = InputFile(Name,2,AeroFlag,GravFlag,1000,Nstep,Nvtk,0,VecNul,1,VecNul,1,self.Wing,Vinf,Rho,AlphaAC,BetaAC,0.,1.)       
        self.Input.AppendTimeFunction(TF)
        self.Input.WriteInputFile()
        self.Input.WriteInitFile()
        # adjustment of the time simulation to catch a flutter sequence
        while (flutter_flag == 1) :
            Nstep = int(Time/DeltaT)
            if verbosity == 1:
                print("simu_time =",Time," s, nstep= ",Nstep)
            # call of gebtaero, the velocity of the input file is override by v= of the command
            Command = "gebtaero -s "+self.Input.GetFileName()+" v="+str(Vinf)+" aero="+str(AeroFlag)+" time="+str(Time)+" vz="+str(Vz)+" nstep="+str(Nstep)+" tfe="+str(PeriodFlutter)+" tfper="+str(PeriodFlutter)+" flutter_limit="+str(FlutterLimit)
            Output = sp.getoutput(Command)
            Output = Output.split()
            if (Output[0] == "#"):
                step_flutter = float(Output[1])
                Time = 0.95*Time*step_flutter/Nstep
            else:
                flutter_flag = 0
            # ~ if (Output[0] == "no"):
                # ~ flutter_flag = 0
            # ~ else:
                # ~ print(Output)
                # ~ step_flutter = float(Output[0])
                # ~ Time = 0.95*Time*step_flutter/Nstep
        self.Input.RemoveInputFile()  
        
        
    def TemporalDynamic(self,Rho,Vinf,AeroFlag,AlphaAC,BetaAC,Time,Nstep,FreqPerturb,TimePerturb,CoefPerturb,CoefVinf,Nvtk,GravFlag=0,verbosity=0):  
        # Computation of vtk files corresponding to a temporal flutter
        import shutil 
        #initial perturbation parameters
        VecNul = np.zeros([3])
        #~ Vz = round(CoefPerturb*Vinf,4)
        Vz = CoefPerturb*Vinf
        #time simulation
        flutter_flag = 1
        # Time Function.
        TF = TimeFunction(1,0.,TimePerturb)
        TF.AppendFunctionEntrieHarmonic(1.,FreqPerturb,0.)
        # Name creation
        Name = self.GetWing().GetName()+"TempDyna"
        # suppression of previous vtk folder
        if (Nvtk > 0):
            vtk_folder = Name+".dat_vtk"
            shutil.rmtree(vtk_folder,True)
        # Input file creation
        self.Input = InputFile(Name,2,AeroFlag,GravFlag,1000,Nstep,Nvtk,0,VecNul,1,VecNul,1,self.Wing,Vinf,Rho,AlphaAC,BetaAC,0.,1.)       
        self.Input.AppendTimeFunction(TF)
        self.Input.WriteInputFile()
        self.Input.WriteInitFile()
        # adjustment of the time simulation to catch a flutter sequence
        if verbosity == 1:
            print("simu_time =",Time," s, nstep= ",Nstep)
        # call of gebtaero, the velocity of the input file is override by v= of the command
        Command = "gebtaero -s "+self.Input.GetFileName()+" v="+str(Vinf)+" aero="+str(AeroFlag)+" time="+str(Time)+" vz="+str(Vz)+" nstep="+str(Nstep)+" tfe="+str(FreqPerturb)+" tfper="+str(TimePerturb)
        Output = sp.getoutput(Command)
        self.Input.RemoveInputFile()  
                
        
    def TemporalDivergenceSpeed(self,Rho,Vmin,Vmax,Vstep,DeltaV,AeroFlag,BetaAC,GravFlag=0,verbosity=0):
        Name = self.GetWing().GetName()+"TemporalDivergence"
        VecNul = np.zeros([3])
        self.Input = InputFile(Name,1,AeroFlag,GravFlag,1,1,0,0,VecNul,0,VecNul,0,self.Wing,Vmin,Rho,1e-6,BetaAC,0.,0.)       
        self.Input.WriteInputFile()
        FileName = self.Input.GetFileName()
        Velocity = Vmin
        
        while Velocity < Vmax and Vstep > DeltaV:
            v = "v="+str(Velocity)
            Command = ["gebtaero","-p",FileName,v]
            Output = ReadLoadsFromPipe(Command)
            Fz = Output[0][2]
            if Velocity == Vmin:
                if Fz > 0.:
                    raise RuntimeError("divergence speed is before Vmin")
                    RemoveFiles()
            if Fz > 0.:                     
                Vstep = 0.5*Vstep
                Velocity = round(Velocity-Vstep,6)
            else:                
                Velocity = round(Velocity+Vstep,6)      
        RemoveFiles()

        if Velocity >= Vmax:
            if verbosity==1:
                print("no divergence speed detectd before Vmax")    
            return -1    
        else:
            if verbosity==1:
                print("divergence speed detected between {0} and {1} m/s".format(round(Velocity-DeltaV,6),Velocity))
            return (Velocity-0.5*DeltaV)        
                    

