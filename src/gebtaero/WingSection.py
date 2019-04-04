# coding=UTF-8

## This class represents a wing section corresponding to a beam element between two keypoints in the fortran solver
class WingSection:
    def __init__(self,Chord,ParameterA,NumberOfElements,SectionLength,CrossSection,Frame):
        ## Chord the wing section chord 
        self.Chord = Chord
        ## The distance betwwen the elstic center and the midchord (in midchord %). Positif if the elastic center is on the side of the trailing edge.
        self.ParameterA = ParameterA
        ## Number of finite element to generate in the wing section
        self.NumberOfElements = NumberOfElements
        ## The length of the wing section
        self.SectionLength = SectionLength
        ## the 
        self.HalfChord = 0.5*Chord
        self.CrossSection = CrossSection
        self.Frame = Frame

    def SetParameterA(self, ParameterA):
            self.ParameterA = ParameterA

    def GetParameterA(self):
            return self.ParameterA

    def SetChord(self, Chord):
        self.Chord = Chord
        self.HalfChord = 0.5*Chord

    def GetChord(self):
        return self.Chord
        
    def GetHalfChord(self):
        return self.HalfChord    


    def SetNumberOfElements(self, NumberOfElements):
        self.NumberOfElements = NumberOfElements
        
    def GetNumberOfElements(self):
        return self.NumberOfElements
        
        
    def SetSectionLength(self,SectionLength):
        self.SectionLength = SectionLength
        
    def GetSectionLength(self):
        return self.SectionLength
        
    def GetCrossSection(self):
        return self.CrossSection
        
    def GetFrame(self):
        return self.Frame
