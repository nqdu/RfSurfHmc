# RFSurfHMC
Joint inversion of Receiver Function and Surface Wave Disperion by Hamiltonian
Monte Carlo Method

Requirements
-------------------
RFSurfHMC has been built upon a few modern packages for its performance and sustainability, as listed below:

name | version|  
|:---:|:----:|  
|[Pybind11](https://github.com/pybind/pybind11) | >=0.4| 
|[FFTW3](https://www.fftw.org/) | >=3.3 |
|[GCC](https://gcc.gnu.org/)| >=7.5|
|[CMAKE](https://cmake.org/)| >=3.10.0|
|mpi4py| >=3.0 |

Installation 
-------------------
+ This package requires some py-package as follow:
    - mpi4py
    - pybind11
    - numpy
    - h5py
    - pyyaml
we recommand install those packages by anaconda
```shell
conda create -n hmc_inv python=3.9
conda install mpi4py pybind11-global numpy h5py pyyaml
```

+ Then use the following to compile the code
```
mkdir -p build
cd build
cmake .. 
make -j4 
make install
```

Usage 
--------------------
+ set your parameters in `param.yaml`
+ For naive HMC, try: `mpiexec -n 4 python main_base.py ` for parallel running
+ For HMC with dual averaging ([hoffman 2014](https://jmlr.org/papers/volume15/hoffman14a/hoffman14a.pdf),algorithm 5), try `mpiexec -n 4 python main_DA.py `
+ `python plot.py` for drawing figures.

New Features
-------------------
+ Python extensions for surface wave dispersion, receiver functions and their Frechet kernels.


+ A more general framework for HMC. 

Update 
-------------------
+ 2023-08-23 we update the forward calculation of receiver function. RF could be calculated in time domain or frequency domain either.
+ 2025-07-26 update pybind11 interface, and add dual averaging, which significantly reduces the need for manual hyperparameter tunings on `dt` and `L`.

References 
-------------------
Junliu Suwen, Qi‐Fu Chen, Nanqiao Du; Joint Inversion of Receiver Function and Surface Wave Dispersion by Hamiltonian Monte Carlo Sampling. Seismological Research Letters 2022;; 94 (1): 369–384. doi: https://doi.org/10.1785/0220220044
