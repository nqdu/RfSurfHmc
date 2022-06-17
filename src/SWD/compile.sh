CC="g++ -O3 -march=native -fPIC -shared"
FC="gfortran -O3 -march=native -fPIC -shared "
include=`python -m pybind11 --includes`
FFLAGS="-ffixed-line-length-none -ffloat-store -w -m64 -mcmodel=medium"
outname="../libsurf"`python3-config --extension-suffix`

set -x 
$FC $FFLAGS -g -c surfdisp96.f -o surfdisp96.o
$CC -g -c surfdisp.cpp -o surfdisp.o 
$FC -g -c slegn96.f90 -o slegn96.o 
$FC -g -c sregn96.f90 -o sregn96.o 
$CC -g -c main.cpp -o main.o  $include
$CC main.o surfdisp.o sregn96.o slegn96.o surfdisp96.o -o ../$outname -lgfortran
rm *.o *.mod