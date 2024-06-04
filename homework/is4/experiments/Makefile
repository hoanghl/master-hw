CC = g++
FLAGS =-march=native -std=c++2a -mavx -fopenmp -O3
#-fopenmp -O3
BIN = bin

e1: experiment1.cc exp
	${CC} ${FLAGS} experiment1.cc exp.o -o ${BIN}/e1 

exp: exp.cc
	${CC} ${FLAGS} -c exp.cc

clean:
