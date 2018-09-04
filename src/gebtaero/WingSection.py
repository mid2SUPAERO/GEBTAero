# coding=UTF-8

class WingSection:
    """
    Section of a wing. A section is defined for each changing parameter (dihedral,
    flexibility,...)

    :version: 12/02/18
    :author: Bertrand Kirsch
    
    ATTRIBUTES

    The loccal chord of the wing section (m)
    Chord  (private) : float

    Distance between the moment calculation point and the half chord. For isotropic
    wing the moment calculation point is placed on the elastic center
    ParameterA  (private) : float

    Half of the local wing chord used in the aerodynamic model
    HalfChord  (private) : float

    Number of FE elements of the wing section
    NumberOfElements  (private) : int

    Length of the wing section
    SectionLength (private) : float

    """
    def __init__(self,Chord,ParameterA,NumberOfElements,SectionLength,CrossSection,Frame):
        self.Chord = Chord
        self.ParameterA = ParameterA
        self.NumberOfElements = NumberOfElements
        self.SectionLength = SectionLength
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
