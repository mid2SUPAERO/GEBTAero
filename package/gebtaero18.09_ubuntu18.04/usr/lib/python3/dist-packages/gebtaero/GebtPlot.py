# coding=UTF-8
import matplotlib.pyplot as plt
import numpy as np
import os as os
from .utils import *

class GebtPlot:

    """
    Class containing ploting routines
    """
    def EigenFreqDamping(Velocity,Modes,Style=None,ReducedDamping=True,DampAxis=False,DampScale=0.01):
        """
        Plot frequencies and Damping of a set of computed modes
        """
        Nmodes = len(Modes[0,:,0])
        plt.figure(1)
        plt.subplot(211)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.ylabel('frequencies (Hz)',fontsize=30)
        for i in range(Nmodes):
            if Style is None:
                plt.plot(Velocity,Modes[:,i,1],linewidth=3)
            else:    
                plt.plot(Velocity,Modes[:,i,1],Style,linewidth=3)
        # ~ plt.axis([min(Velocity),max(Velocity),0.,1.1*np.amax(Modes[:,:,1])])
        # ~ deg5, = plt.plot(Velocity,Modes[:,0,1],linewidth=3,label="[-45/2$z_c$/5]")
        # ~ deg10, = plt.plot(Velocity,Modes[:,1,1],linewidth=3,label="[-45/2$z_c$/10]")
        # ~ plt.axis([min(Velocity),max(Velocity),25.,45.])
        # ~ plt.legend(handles=[deg5,deg10],loc=1, prop={'size': 25})
        plt.grid()		

        plt.subplot(212)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlabel('aerodynamic speed (m/s)',fontsize=30)
        if ReducedDamping:
            plt.ylabel('reduced damping',fontsize=30)
            for i in range(Nmodes):
                for j in range(len(Modes[:,i,0])):
                    if Modes[j,i,1] is not None:
                        if abs(Modes[j,i,1]) <1e-12 :
                            Modes[j,i,0] = None
                        else :
                            Modes[j,i,0] = Modes[j,i,0] /(-4*np.pi*Modes[j,i,1])
                if Style is None:
                    plt.plot(Velocity,Modes[:,i,0],linewidth=3)
                else:    
                    plt.plot(Velocity,Modes[:,i,0],Style,linewidth=3)
        else :
            plt.ylabel('real part',fontsize=30)
            for i in range(Nmodes):
                if Style is None:
                    plt.plot(Velocity,Modes[:,i,0],linewidth=3)
                else:    
                    plt.plot(Velocity,Modes[:,i,0],Style,linewidth=3)            
        if DampAxis:    
            plt.axis([min(Velocity),max(Velocity),-DampScale,DampScale])    
        plt.hlines(0.0, min(Velocity),max(Velocity), colors='r', linestyles='dashdot', label='', hold=None, data=None,linewidth=3)
        plt.grid()
        plt.show()
        
        
    def EigenFreqDampingUnsorted(Velocity,Modes,Style=None,DampAxis=False,DampScale=0.01):
        """
        Plot frequencies and Damping of a set of computed modes
        """
        plt.figure(1)
        plt.subplot(211)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.ylabel('frequencies (Hz)',fontsize=30)
        for i in range(len(Velocity)):
            temp = Modes[i]
            for j in range(len(temp)):
                plt.plot(Velocity[i],temp[j,1],"k.",linewidth=2)
        plt.grid()		

        plt.subplot(212)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlabel('aerodynamic speed (m/s)',fontsize=30)
        plt.ylabel('reduced damping',fontsize=30)
        for i in range(len(Velocity)):
            temp = Modes[i]
            for j in range(len(temp)):
                plt.plot(Velocity[i],temp[j,0],"k.",linewidth=2)
        if DampAxis:    
            plt.axis([min(Velocity),max(Velocity),-DampScale,DampScale])    
        plt.grid()
        plt.show()

    def WriteParaviewScript(VtkFolder,VtkName,AirfoilPath="/opt/gebtaero/airfoil/airfoil_default.vtk"):
        file = open('pvscript.py',"w")
        NbVtk = len(os.listdir(VtkFolder))
        FileList = []
        for i in range(NbVtk):
            FileList.append(VtkFolder+VtkName+str(i+1)+'.vtk')
        
        file.write('from paraview.simple import *\n')
        file.write('AirfoilVtkReader = LegacyVTKReader(FileNames=["'+str(AirfoilPath)+'"])\n')
        file.write("renderView1 = GetActiveViewOrCreate('RenderView')\n")
        file.write("transform1 = Transform(Input=AirfoilVtkReader)\n")
        file.write("transform1.Transform.Rotate = [90.0, 0.0, 0.0]\n")
        file.write("WingVtkReader = LegacyVTKReader(FileNames="+str(FileList)+")\n")
        file.write("animationScene1 = GetAnimationScene()\n")
        file.write("animationScene1.UpdateAnimationUsingDataTimeSteps()\n")
        file.write("warpByVector1 = WarpByVector(Input=WingVtkReader)\n")
        file.write("warpByVector1.Vectors = ['POINTS', 'displacements']\n")
        file.write("glyphWithCustomSource1 = GlyphWithCustomSource(Input=warpByVector1,GlyphType=transform1)\n")
        file.write("glyphWithCustomSource1.Scalars = [None, '']\n")
        file.write("glyphWithCustomSource1.Vectors = ['POINTS', 'midchord_ori']\n")
        file.write("glyphWithCustomSource1.ScaleFactor = 1.\n")
        file.write("glyphWithCustomSource1.ScaleMode = 'vector'\n")
        file.write("glyphWithCustomSource1Display = Show(glyphWithCustomSource1, renderView1)\n")
        file.write("glyphWithCustomSource1Display.DiffuseColor = [0.5, 0.0, 0.0]\n")
        file.write("renderView1.Update()\n")
        file.write("renderView1.InteractionMode = '2D'\n")
        file.write("renderView1.CameraPosition = [10,-15,-5]\n")
        file.write("renderView1.CameraFocalPoint = [0.0, 0.0,0.0]\n")
        file.write("renderView1.CameraViewUp = [0,0,-1]\n")
        file.write("renderView1.ResetCamera()\n")
        
    def ParaviewOutput(self,VtkFolder,VtkName,AirfoilPath="/opt/gebtaero/airfoil/airfoil_default.vtk"):
        self.WriteParaviewScript(AirfoilPath,VtkFolder,VtkName)
        RunParaviewScript(pvscript.py)
        
