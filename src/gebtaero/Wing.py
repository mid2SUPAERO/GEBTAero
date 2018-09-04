# coding=UTF-8
import numpy as np

class Wing:
    """
    Wing compose by one ore more wing section
    :version: 12/02/18
    :author: Bertrand Kirsch
    
    ATTRIBUTES
    
    A list containing the wing sections
    WingSections (private) : WingSection[]
    
    A vector containing the wing root position
    WingRootPosition (private) : np.array[3]
    
    A vector containing the direction of the wing at the root
    WingRootAxis (private) : np.array[3]
    
    The wing twist at the root (rad)
    WingRootTwist (private) : float
    
    The wing frame to be input in the solver regarding the wing axis and the wing twist
    WingFrame (private) : np.array[3][3]        
    """
    def __init__(self,Name,WingRootPosition):
        self.Name = Name
        self.WingRootPosition = np.copy(WingRootPosition)
        
        # Initialisation of list
        self.WingSections = []
        self.Frames = []
        self.CrossSections = []
        self.KpList = []
        self.KpList.append(WingRootPosition)
        
    def GetWingRootPosition(self):
        return self.WingRootPosition
    
    def GetName(self):
        return self.Name    
        
    def AppendWingSection(self,WingSection):
        """
        Add a section to the current wing, first section is at the wing root, last section at the wing tip.
        """
        # add the wing section
        self.WingSections.append(WingSection)
        
        # add the cross section if not in the list
        if WingSection.GetCrossSection() not in self.CrossSections:
            self.CrossSections.append(WingSection.GetCrossSection())
        # add the frame if not in the list 
        if WingSection.GetFrame() not in self.Frames:
            self.Frames.append(WingSection.GetFrame())
        # add a keypoint to the list
        self.KpList.append(self.KpList[-1]+WingSection.GetFrame().GetAxis()*WingSection.GetSectionLength())
           
    def GetWingSections(self):
        return self.WingSections
        
    def GetCrossSections(self):
        return self.CrossSections

    def GetFrames(self):
        return self.Frames
        
    def GetKpList(self):
        return self.KpList    
        
    def GetWeight(self):
        WS = self.GetWingSections()
        Weight = 0.
        for Section in WS:
            Weight = Weight + Section.GetSectionLength()*Section.GetCrossSection().GetMassMatrix()[0][0]
        Weight = Weight*9.81
        return Weight     
        
    def GetSurface(self):
        WS = self.GetWingSections()
        Surface = 0.
        for Section in WS:
            Surface = Surface + Section.GetSectionLength()*Section.GetChord()
        return Surface    
