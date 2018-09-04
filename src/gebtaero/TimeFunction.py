# coding=UTF-8

class TimeFunction:
    """
    Time function object as defined in gebt documentation
    """
    def __init__(self,FunctionType,FunctionStart,FunctionEnd):
        self.FunctionType = FunctionType
        self.FunctionStart = FunctionStart
        self.FunctionEnd = FunctionEnd
        
        # Initialisation of function entries list
        self.FunctionEntries = []
    
    def GetFunctionType(self):
        return self.FunctionType
        
    def GetFunctionStart(self):
        return self.FunctionStart
        
    def GetFunctionEnd(self):
        return self.FunctionEnd
        
    def AppendFunctionEntriePieceWise(self,time,value):
        if(self.FunctionType == 0):
            self.FunctionEntries.append([time,value])           
        else:
            raise TypeError('the time function is not piecewise')
    
    def AppendFunctionEntrieHarmonic(self,amplitude,period,phase):
        if(self.FunctionType == 1):
            self.FunctionEntries.append([amplitude,period,phase])           
        else:
            raise TypeError('the time function is not piecewise')
            
    def GetFunctionEntrie(self,index):
        return self.FunctionEntries[index]
    
    def GetFunctionEntries(self):
        return self.FunctionEntries
