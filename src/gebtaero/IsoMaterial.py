# coding=UTF-8

## This class correspond to the TYPE=ISO of the *ELASTIC input card of ccx solver
class IsoMaterial:
    def __init__(self,E,Nu,Rho=None):
        ## Young Modulus
        self.E = E
        ## Poisson Ratio
        self.Nu = Nu
        ## Density
        self.Rho = Rho
        
    def GetIso(self):
        return [self.E,self.Nu]    
    
    def GetDensity(self):
        return self.Rho    
