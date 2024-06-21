CC = clang++
BIN = bin
FLAGS = -march=native -std=c++2a -fopenmp -O3
# -O3

e1: exp1.cc exp.o
	${CC} ${FLAGS} exp1.cc exp.o -o ${BIN}/e1

exp.o: exp.cc
	${CC} ${FLAGS} exp.cc