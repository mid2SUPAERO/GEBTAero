class OrthoMaterial:
    """
    This class contains the orthotropic material characteristics
    """
    def __init__(self,El,Et,Nult,Glt,Rho=None):
        self.El = El
        self.Et = Et
        self.Nult = Nult
        self.Glt = Glt
        self.Rho = Rho
        
    def GetOrtho(self):
        return [self.El,self.Et,self.Nult,self.Glt]    
    
    def GetDensity(self):
        return self.Rho    
