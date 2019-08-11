# GEBTAero 

GEBTAero is an aeroelasticity simulation toolbox with a computation code in Fortran and a pre/postprocessor coded in Python.
It is dedicated to quick critical speed computation, particularly for aeroelastic tailoring optimisation of high aspect ratio composite wing.
The computation code is derived from GEBT program developped by Prof. Yu (https://cdmhub.org/resources/gebt).
The pre/postprocessor uses several open source programs available in most linux distros repositories:
* Calculix : a Finite Element Method solver (http://www.calculix.de/)
* Paraview : a  data analysis and visualization application (https://www.paraview.org/)
* Mumps : a parallel sparse direct solver (http://mumps.enseeiht.fr/)
* Arpack : a sparse eigenvalue solver (https://www.caam.rice.edu/software/ARPACK/)

## Installation
Two options are available:

### Debian package
For Ubuntu 18.04, 19.04 and Debian 10, download the .deb file available in the package folder and launch it.
It will automatically install all the dependancies.
For other linux distributions, you can ask for a .deb or .rpm package creation.

### Compilation

Install the dependancies. 
On Ubuntu:
```bash
sudo apt install paraview calculix-ccx calculix-cgx libmumps-seq-dev libarpack2-dev python3 python3-numpy python3-matplotlib gfortran make git doxygen
```

Clone the repository ("/MyFolderPath/" is the path to the folder where the  GEBTAero folder will be copied)
```bash
cd /MyFolderPath/
git clone https://framagit.org/BertrandK/GEBTAero.git
```

Compile gebtaero (fortran executable of the computation core) and unical (mesh format translator from unv to inp).
You can skip that step if you prefer to use already compile executable.
```bash
cd /MyFolderPath/bin/
make clean
make gebtaero
make unical
```

Create a symbolink link between to the python script folder and the two executable (gebtaero and unical)
For example, with Ubuntu : 
```bash
cd /usr/lib/python3/dist-packages/
sudo ln -s /MyFolderPath/GEBTAero/src/gebtaero/
cd /usr/bin/
sudo ln -s /MyFolderPath/GEBTAero/bin/gebtaero
sudo ln -s /MyFolderPath/GEBTAero/bin/unical
```

## Testing

The folder cas_test is a set of automated python script designed to test many program functionalities.
in a terminal launch:
```bash
cd /MyFolderPath/GEBTAero/cas_test/
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
It will generate a .ech text file summarising the input parameters, a .out text file with the output data and optionally a vtk file folder

You can use an IDE to easily launch the scripts and read the source code. 
This program has been developed using geany (www.geany.org). Open a python script and press "Maj+F5". Modify "python" in "python3" and press "F5" to launch the script.

## Using GEBTAero within OpenDMAO framework
GEBTAero solver is designed to be used within a multidisciplinary optimisation framework like openMDAO (https://openmdao.org/). To illustrate that point, a simple optimisation application is proposed in examples/optimisation/ folder.
To install openMDAO on Ubuntu or Debian :
```bash
sudo apt install python3-pip
pip3 install openmdao
```

## Windows users
GEBTAero has been developped and tested on linux OS. If you are a Windows user :

### Windows compilation
You can modify the program to compile it on Windows. The main difficulties concern external programs and libraries used like Calculix and Mumps.

### Virtual machine
Using VirtualBox (https://www.virtualbox.org/) or Vmware Workstation Player (https://www.vmware.com/) and a Ubuntu 18.04 iso (https://ubuntu.com/download/desktop), you can use GEBTAero as a linux user.

### WSL
Using WSL (Windows Subsystem for Linux) on a compatible Windows 10 version (https://docs.microsoft.com/windows/wsl/install-win10), you can install within the Windows shop a Ubuntu 18.04 kernel and install GEBTAero on it.
Then you just need to had "wsl" in windows bash command before linux command.

## Documentation
If you have installed the debian package, an application "GEBTAero Doc" is installed on your system. It launches an html page with the whole documentation of the sources (Fortran and Python)
You can also generate the documentation using doxygen and the file "Doxyfile" in the "doc" folder.
```bash
cd /MyFolderPath/GEBTAero/doc/
doxygen Doxyfile
```
To launch the documentation main page :
```bash
cd /MyFolderPath/GEBTAero/doc/html/
firefox index.html
```

## License
See the License file in the repository

## Acknowledgement
This research work was funded by the French Air Force Academy Research Center in collaboration with ISAE-Supaero



