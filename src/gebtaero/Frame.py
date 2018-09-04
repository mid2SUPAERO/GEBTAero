# coding=UTF-8
import numpy as np
import math

class Frame:
    """
    class containing the matrix Cab of the solver ie the local frame of the undeformed wing
    """    
    def __init__(self,Axis,Twist):
        self.Axis = Axis
        self.Twist = Twist
        self.FrameMatrix = np.zeros([3,3])
        # determination of the frame (frame b in the solver)
        # Frame a
        Xa = np.zeros([3])
        Xa[0] = 1.
        Ya = np.zeros([3])
        Ya[1] = 1.
        Za = np.zeros([3])
        Za[2] = 1.
        # ~ # Xb = (WingAxis.Ya).Ya + (WingAxis.Za).Za
        # ~ if(self.Axis[1]>0):
            # ~ Xb= np.vdot(self.Axis,Ya)*Ya+np.vdot(self.Axis,Za)*Za
        # ~ else:
            # ~ Xb= -np.vdot(self.Axis,Ya)*Ya-np.vdot(self.Axis,Za)*Za
        # ~ Xb = Xb/np.sqrt(np.vdot(Xb,Xb))
        # Xb = +/-WingAxis
        if(self.Axis[1]>0):
            Xb= Axis

        else:
            Xb= -Axis
        Yb = np.cross(Xb,Za)
        Yb = 1./np.linalg.norm(Yb)*Yb
        Zb = np.cross(Xb,Yb)
        # Euler-Rodrigues formula (cf wikipedia) to make the rotation of WingTwist about the Wing axis

        a = math.cos(0.5*self.Twist)
        w = Xb*math.sin(0.5*self.Twist)

        # Yb is the rotation of Yb with an angle WingTwist about Xb 
        Yb  = Yb + 2*a*np.cross(w,Yb)+2*np.cross(w,np.cross(w,Yb))
        # Zb is the rotation of Zb with an angle WingTwist about Xb
        Zb = Zb + 2*a*np.cross(w,Zb)+2*np.cross(w,np.cross(w,Zb))
        #temp
        
        
        # ~ self.FrameMatrix[0][:] = Xb
        # ~ self.FrameMatrix[1][:] = Yb
        # ~ self.FrameMatrix[2][:] = Zb

        self.FrameMatrix[:,0] = Xb
        self.FrameMatrix[:,1] = Yb
        self.FrameMatrix[:,2] = Zb

    def GetFrameMatrix(self):
        return self.FrameMatrix
        
    def GetAxis(self):
        return self.Axis
        
    def GetTwist(self):
        return self.Twist         
