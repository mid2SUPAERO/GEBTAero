# GEBTAero 

GEBTAero is an aeroelasticity simulation toolbox with a computation code coded in Fortran and a pre/postprocessor coded in Python.
The computation code is derived from GEBT program developped by Prof. Yu (https://cdmhub.org/resources/gebt).
The pre/postprocessor uses several open source programs available in most linux distros repositories:
* calculix : a Finite Element Method solver (http://www.calculix.de/)
* paraview : a  data analysis and visualization application (https://www.paraview.org/)
* Mumps : a parallel sparse direct solver (http://mumps.enseeiht.fr/)
* Arpack : a sparse eigenvalue solver (https://www.caam.rice.edu/software/ARPACK/)

## Installation
Two options are available:

### Debian package
For Ubuntu 18.04 and Debian 10, download the .deb file and launch it.
It will automatically install all the dependancies.
For other linux distributions, you can ask for a .deb or .rpm package creation.

### Compilation
Clone the repository 
use the MakeFile in the bin folder and adapt it to your system.

Install the dependancies. On Ubuntu:
```bash
sudo apt install paraview calculix-ccx calculix-cgx libmumps-seq-dev libarpack2-dev python3 python3-numpy python3-matplotlib
```
Compile gebtaero and unical (mesh format translator from unv to inp)


## Testing

The folder cas_test is a set of automated python script designed to test many program functionalities.
in a terminal launch:
```bash
python3 tests.py
```

## Usage
Besides cas_test folder, examples folder contains a set of detailled script designed to help you to set your own problems.
The pre/postprocessor script must be launch with python3 (not python2).
```bash
python3 myscript.py
```
You can also directly use the computation code with .dat file (show examples):
```bash
gebtaero example.dat
```

## Licence
GNU General Public License v3.0



