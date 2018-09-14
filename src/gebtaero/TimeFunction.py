# coding=UTF-8

## This class interface the concept of time function timefunctionmodule::timefunction with the fortran solver
class TimeFunction:
    def __init__(self,FunctionType,FunctionStart,FunctionEnd):
        ## link to timefunctionmodule::timefunction::fun_type
        self.FunctionType = FunctionType
        ## link to timefunctionmodule::timefunction::ts
        self.FunctionStart = FunctionStart
        ## link to timefunctionmodule::timefunction::te        
        self.FunctionEnd = FunctionEnd
        
        # Initialisation of function entries list
        self.FunctionEntries = []
    
    def GetFunctionType(self):
        return self.FunctionType
        
    def GetFunctionStart(self):
        return self.FunctionStart
        
    def GetFunctionEnd(self):
        return self.FunctionEnd
        
    ##Append a function entrie to a piecewise function
    #@param time link to timefunctionmodule::timefunction::time_val
    #@param value link to timefunctionmodule::timefunction::fun_val
    def AppendFunctionEntriePieceWise(self,time,value):
        if(self.FunctionType == 0):
            self.FunctionEntries.append([time,value])           
        else:
            raise TypeError('the time function is not piecewise')
    
    ##Append a function entrie to an harmonic function
    #@param link to timefunctionmodule::timefunction::time_val
    #@param link to timefunctionmodule::timefunction::fun_val
    #@param link to timefunctionmodule::timefunction::phase_val
    def AppendFunctionEntrieHarmonic(self,amplitude,period,phase):
        if(self.FunctionType == 1):
            self.FunctionEntries.append([amplitude,period,phase])           
        else:
            raise TypeError('the time function is not piecewise')
            
    def GetFunctionEntrie(self,index):
        return self.FunctionEntries[index]
    
    def GetFunctionEntries(self):
        return self.FunctionEntries
