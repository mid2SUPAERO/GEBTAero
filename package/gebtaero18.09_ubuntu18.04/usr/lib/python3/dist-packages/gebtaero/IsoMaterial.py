class IsoMaterial:
    """
    This class contains the isotropic material characteristics
    """
    def __init__(self,E,Nu,Rho=None):
        self.E = E
        self.Nu = Nu
        self.Rho = Rho
        
    def GetIso(self):
        return [self.E,self.Nu]    
    
    def GetDensity(self):
        return self.Rho    
