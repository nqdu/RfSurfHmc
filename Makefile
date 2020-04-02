cc = g++
cc1 = gfortran
prom1 = bin/synthetic
prom2 = bin/inversion
opt=-O3 -std=c++11

f77flags= -ffixed-line-length-none -ffloat-store -W  -fbounds-check
deps = $(shell find include/*.h)
cppsrc = $(shell find src/*.cc)
f77src = $(shell find src/*.f)
objcpp = $(cppsrc:%.cc=%.o) 
objf77 =$(f77src:%.f=%.o)
synobj=src/surfdisp96.o src/disp.o src/gaussian.o src/synthetic.o 
invobj=src/surfdisp96.o src/celldisp.o src/disp.o src/gaussian.o src/main.o

all: $(prom1) $(prom2)
$(prom1):$(synobj)
	$(cc1) -o $(prom1) $(synobj) $(opt) -lstdc++ 
$(prom2):$(invobj)
	$(cc1) -o $(prom2) $(invobj) $(opt) -lstdc++

%.o: %.cc
	$(cc) -g -c -Iinclude/ $< -o $@ $(opt) -lm
%.o:%.f
	$(cc1) -g -c $< -o $@ $(f77flags) -O3
clean:
	rm $(objcpp) $(objf77)

