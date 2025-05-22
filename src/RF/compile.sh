CC="g++ -O3 -march=native -fPIC -shared"
FC="gfortran -O3 -march=native -fPIC -shared"
#FC="gfortran  -fPIC"
include=`python -m pybind11 --includes`

FFTW_INC="-I/opt/fftw-3.3.10_share/include"
FFTW_LIB="-L/opt/fftw-3.3.10_share/lib -lfftw3"
outname="../../librf"`python3-config --extension-suffix`


set -x 
$FC -g -c RFModule.f90 -o RFModule.o $FFTW_INC 
$FC -g -c deconit.f90 -o deconit.o $FFTW_INC 
$FC -g -c fftpack.f90 -o fftpack.o $FFTW_INC 
$CC -g -c main.cpp -o main.o $include
$CC main.o RFModule.o deconit.o fftpack.o -o ./$outname $FFTW_LIB -lgfortran

rm *.o *.mod

# there might be some bugs for -march=native 
