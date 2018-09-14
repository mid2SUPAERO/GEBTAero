# coding=UTF-8

## This class define a composite ply 
class CompositePly:

    def __init__(self,Material,Thickness,Orientation):
        ## the material of the ply (either IsoMaterial or OrthoMaterial)
        self.Material = Material
        ## the thickness of the ply (m)
        self.Thickness = Thickness
        ## the orientation of fiber in case of OrthoMaterial
        self.Orientation = Orientation
        
    def GetMaterial(self):
        return self.Material
        
    def GetThickness(self):
        return self.Thickness
        
    def GetOrientation(self):
        return self.Orientation
