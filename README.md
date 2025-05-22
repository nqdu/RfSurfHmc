# RFSurfHMC
Joint inversion of Receiver Function and Surface Wave Disperion by Hamiltonian
Monte Carlo Method

# Pre-update for v1.0
1. NUTS (No-U-Turn Sampling) will be implemented
2. Add taper for receiver function.
3. code reorganize

Requirements
-------------------
RFSurfHMC has been built upon a few modern packages for its performance and sustainability, as listed below:

name | version|  
|:---:|:----:|  
|[Pybind11](https://github.com/pybind/pybind11) | >=0.4| 
|[FFTW3](https://www.fftw.org/) | >=3.3 |
|[GCC](https://gcc.gnu.org/)| >=7.5|
|mpi4py| >=3.0 |

Installation 
-------------------
+ This package requires some py-package as follow:
    - mpi4py
    - pybind11
    - numpy
we recommand install those packages by anaconda
```shell
conda create -n hmc_inv python=3.7
conda install mpi4py pybind11 numpy 
```

+ This package requires fftw also. Download and install the fftw by the instruction fftw.org
Modify the compile.sh in `./src/RF` and `./src/SWD`, and subsitude the var `FFTW_INC` and `FFTW_LIB`

+ Then run `compile.sh` to run this code

Usage 
--------------------

+ execute `mpiexec -n 16 python main_joint.py ` for parallel running.

or `python main_joint.py` for serial running.

+ `python plot.py` for drawing figures.

New Features
-------------------
+ Python extensions for surface wave dispersion, receiver functions and their Frechet kernels.
+ A more general framework for HMC. 

References 
-------------------
TBD
