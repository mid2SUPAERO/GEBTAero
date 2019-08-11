from openmdao.api import Problem,IndepVarComp,ScipyOptimizeDriver
from Classes   import SouplesseMat
from time import *
import os
name = ctime()
os.makedirs(name)
os.chdir(name)

objectif='mCoupling' # opposite of the bend/twist coupling coefficient
objectif='coupling'  # bend/twist coupling coefficient
objectif='bendflex' # bending flexibility
objectif='twistflex' # twisiting flexibility

Start = time()
Prob = Problem()
model = Prob.model
model.add_subsystem('Souplessesyst',SouplesseMat.SouplesseMat(),promotes=['*'])
model.add_subsystem('IndepVar1',IndepVarComp('Ori_1',0,units='deg'),promotes=['Ori_1'])
model.add_subsystem('IndepVar2',IndepVarComp('Ori_2',0,units='deg'),promotes=['Ori_2'])


## run optimisation

## COBYLA
Prob.driver = ScipyOptimizeDriver() 
Prob.driver.options['optimizer'] = 'COBYLA' 

Prob.driver.options['maxiter'] = 1000
Prob.driver.options['tol'] = 1e-9

model.add_design_var('Ori_1', lower=-90, upper=90)
model.add_design_var('Ori_2', lower=-90, upper=90) 
model.add_objective(objectif) 



Prob.setup() 
Prob.run_driver() 

print(Prob['Ori_1'])
print(Prob['Ori_2'])
print(Prob[objectif])

End = time()

print("Elapsed time : {0}".format(End-Start))

#creation of a log file
file = open("logfile.txt","w")
file.write("Optimisation COBYLA\n")
file.write(name+"\n")
file.write("Orientation pli inférieur = {0}\n".format(Prob['Ori_1']))
file.write("Orientation pli supérieur = {0}\n".format(Prob['Ori_2']))
file.write("Fonction objectif ({0}) = {1}\n".format(objectif,Prob[objectif]))

file.write("Elapsed time : {0}".format(End-Start))
