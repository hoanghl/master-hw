CC = clang++
FLAGS = -g -march=native -std=c++2a
#-fopenmp -O3
BIN = bin

e1: exp_1.cc exp
	${CC} ${FLAGS} exp_1.cc exp.o -o ${BIN}/e1 

e2: exp_2.cc
	${CC} ${FLAGS} exp_2.cc -o ${BIN}/e2 

exp: exp.cc
	${CC} ${FLAGS} -c exp.cc

clean:
