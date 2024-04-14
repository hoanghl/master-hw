
CC=clang++
MPICC=mpicc
FORT=gfortran

BIN=bin

FLAGS = -g -O2

serial_cpp: clean
	${CC} ${FLAGS} -o ${BIN}/md2d_ser_cpp md2d.cpp

serial_for: clean
	${FORT} -O2 -o ${BIN}/md2d_ser_for md2d.f90

clean:
	rm -rf bin/*