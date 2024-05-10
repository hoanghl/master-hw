
CC=clang++
MPICC=mpicc
FORT=gfortran

BIN=bin

FLAGS = -std=c++2a -g

serial_cpp: 
	${CC} ${FLAGS} -o ${BIN}/md2d_ser_cpp md2d.cpp

serial_for: 
	${FORT} -O2 -o ${BIN}/md2d_ser_for md2d.f90

md: md.cc
	${CC} md.cc ${FLAGS} -o ${BIN}/md

clean:
	rm -rf bin/*