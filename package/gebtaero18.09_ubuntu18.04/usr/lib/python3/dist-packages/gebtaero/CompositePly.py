class CompositePly:
    """
    This class defined a laminated composite ply with an orthotropic material, 
    a thickness (m) and a fiber orientation (°)
    """
    def __init__(self,Material,Thickness,Orientation):
        self.Material = Material
        self.Thickness = Thickness
        self.Orientation = Orientation
        
    def GetMaterial(self):
        return self.Material
        
    def GetThickness(self):
        return self.Thickness
        
    def GetOrientation(self):
        return self.Orientation
