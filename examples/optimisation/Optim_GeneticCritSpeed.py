from openmdao.api import Problem,IndepVarComp,SimpleGADriver,ParallelGroup
from Classes   import Vflutter
from time import *
import os
name = ctime()
os.makedirs(name)
os.chdir(name)

objectif='mVflutter'
Ori1Dep = -42
Ori2Dep = -42

Start = time()
Prob = Problem()
model = Prob.model
model.add_subsystem('Souplessesyst',Vflutter.Vflutter(),promotes=['*'])
model.add_subsystem('IndepVar1',IndepVarComp('Ori_1',0,units='deg'),promotes=['Ori_1'])
model.add_subsystem('IndepVar2',IndepVarComp('Ori_2',0,units='deg'),promotes=['Ori_2'])


## run optimisation

## SimpleGA
Prob.driver = SimpleGADriver() 


# ~ Prob.driver.options['bits'] = {'Ori_1': 8}
# ~ Prob.driver.options['bits'] = {'Ori_2': 8}
# ~ Prob.driver.options['run_parallel'] = True
# ~ Prob.driver.options['max_gen'] = 10
# ~ Prob.driver.options['tol'] = 1e-9

model.add_design_var('Ori_1', lower=-90, upper=90)
model.add_design_var('Ori_2', lower=-90, upper=90) 
model.add_objective(objectif) 
# ~ model.add_objective('mCoupling') 

Prob.driver.options['run_parallel'] = True
# ~ Prob.driver.options['procs_per_model'] = 6
# ~ Prob.driver.options['bits'] = {'Ori_1': 6}
# ~ Prob.driver.options['bits'] = {'Ori_2': 6}
Prob.driver.options['max_gen'] = 20

Prob.setup(check=False)
Prob.run_driver() 

print(Prob['Ori_1'])
print(Prob['Ori_2'])
print(Prob[objectif])

End = time()

print("Elapsed time : {0}".format(End-Start))

#creation of a log file
file = open("logfile.txt","w")
file.write("Optimisation genetic\n")
file.write(name+"\n")
file.write("Orientation départ 1 = {0}\n".format(Ori1Dep))
file.write("Orientation départ 2 = {0}\n".format(Ori2Dep))
file.write("Orientation pli inférieur = {0}\n".format(Prob['Ori_1']))
file.write("Orientation pli supérieur = {0}\n".format(Prob['Ori_2']))
file.write("Fonction objectif ({0}) = {1}\n".format(objectif,Prob[objectif]))

file.write("Elapsed time : {0}".format(End-Start))



