CC="g++ -O3 -march=native -fPIC -shared"
FC="gfortran -O3 -march=native -fPIC -shared"
#FC="gfortran  -fPIC"
include=`python -m pybind11 --includes`
FFTW_INC="-I/opt/fftw-3.3.10/include"
FFTW_LIB="-L/opt/fftw-3.3.10/lib -lfftw3"
outname="../../librf"`python3-config --extension-suffix`


set -x 
$FC -g -c RFModule.f90 -o RFModule.o $FFTW_INC 
$CC -g -c main.cpp -o main.o $include
$CC main.o RFModule.o -o ./$outname $FFTW_LIB -lgfortran

rm *.o *.mod

# there might be some bugs for -march=native 
