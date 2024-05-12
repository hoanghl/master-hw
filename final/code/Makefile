
CC=clang++
MPICC=mpic++
FORT=gfortran

BIN=bin

FLAGS = -std=c++2a -O3

md2d: 
	${CC} ${FLAGS} md2d.cpp -o ${BIN}/md2d

serial_for: 
	${FORT} -O3 -o ${BIN}/md2d_ser_for md2d.f90

md2d_parallel: md2d_parallel.cc
	${MPICC} ${FLAGS} md2d_parallel.cc -o ${BIN}/md2d_parallel

clean:
	rm -rf bin/*