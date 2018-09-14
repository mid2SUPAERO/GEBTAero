# coding=UTF-8
import numpy as np
import math

## This class describes the frame b of a wing section (see GEBTAero paper)
class Frame:
    def __init__(self,Axis,Twist):
        ## axis xb
        self.Axis = Axis
        ## twist of the wing section
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
        self.FrameMatrix[:,0] = Xb
        self.FrameMatrix[:,1] = Yb
        self.FrameMatrix[:,2] = Zb

    def GetFrameMatrix(self):
        return self.FrameMatrix
        
    def GetAxis(self):
        return self.Axis
        
    def GetTwist(self):
        return self.Twist         
