# coding=UTF-8

## This class correspond to the TYPE=ENGINEERING CONSTANTS of the *ELASTIC input card of ccx solver
class OrthoMaterial:

    def __init__(self,El,Et,Nult,Glt,Rho=None):
        ## Young Modulus in fiber direction
        self.El = El
        ## Young Modulus in transverse direction
        self.Et = Et
        ## Poisson ratio
        self.Nult = Nult
        ## Coulomb Modulus
        self.Glt = Glt
        ## Density
        self.Rho = Rho
        
    def GetOrtho(self):
        return [self.El,self.Et,self.Nult,self.Glt]    
    
    def GetDensity(self):
        return self.Rho    
