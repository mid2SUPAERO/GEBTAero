# coding=UTF-8
import numpy as np
from gebtaero import *
import subprocess as sp
from time import time

# Change TemporalTest to True if you want to run transient dynamic flutter test (Beware, it may take time)
TemporalTest = False
# ~ TemporalTest = True

Start = time()

def RunTestCase(FileName):
    #~ Command = "python3 "+FileName
    Command = ["python3",FileName]
    return sp.run(Command)
    
# Isotropic wing defined using litterature test cases
# Goland Wing
print("Goland Wing")
FileName = "AilesIso/Goland.py"
RunTestCase(FileName)

# Patil Wing
print("Patil Wing")
FileName = "AilesIso/Patil.py"
RunTestCase(FileName)

# periodic beam 3D FEM homogeneisation method
# Disk (analytic)
print("Analytical cross section (disc)")
FileName = "SectionDroite/Disc.py"
RunTestCase(FileName)

# Plate (Flexibility by FEM and Mass analytic)
print("isotropic plate")
FileName = "SectionDroite/Plate.py"
RunTestCase(FileName)

# Box (Flexibility by FEM and Mass analytic)
print("isotropic box")
FileName = "SectionDroite/Box.py"
RunTestCase(FileName)

if TemporalTest:
    # Temporal flutter
    print("Temporal flutter")
    FileName = "AilesIso/Temporal.py"
    RunTestCase(FileName)

print("#####All tests OK#####")

sp.getoutput("rm *.msh")

End = time()
print("------------------------------------------")
print("Elapsed time : {0}".format(End-Start))
