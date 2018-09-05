# GEBTAero ![CI status](https://img.shields.io/badge/build-passing-brightgreen.svg)

GEBTAero is an aeroelasticity simulation toolbox with a computation code coded in Fortran and a pre/postprocessor coded in Python.
The computation code is derived from GEBT developped by Prof. Yu (https://cdmhub.org/resources/gebt).
The pre/postprocessor uses several open source programs available in most linux distros repositories:
* calculix : a Finite Element Method solver (http://www.calculix.de/)
* paraview : a  data analysis and visualization application (https://www.paraview.org/)
* Mumps : a parallel sparse direct solver (http://mumps.enseeiht.fr/)
* Arpack : a sparse eigenvalue solver (https://www.caam.rice.edu/software/ARPACK/)

## Installation
Two options are available
### Debian package
For Ubuntu 18.04 and Debian 10, download the .deb file and launch it.
It will automatically install all the dependancies.

### Compilation
clone the repository 
use the MakeFile in the bin folder


## Testing

The folder cas_test is a set of automated python script designed to test many program functionalities.
in a terminal launch:
```bash
python3 tests.py
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
