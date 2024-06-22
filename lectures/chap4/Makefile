CC = g++
BIN = bin
BUILD = build
SRC = src
FLAGS = -march=native -std=c++2a  -fopenmp -mavx -O3
# -fopenmp -mavx -O3

main: chap4
	${CC} ${FLAGS} -I ./include/ ${BUILD}/* ${SRC}/main.cc -o ${BIN}/main

chap4: 
	${CC} ${FLAGS} -I ./include/ -c ${SRC}/cpu.cc -o ${BUILD}/cpu.o
	${CC} ${FLAGS} -I ./include/ -c ${SRC}/chap4.cc -o ${BUILD}/chap4.o

clean:
	rm -f ${BUILD}/* ${BIN}/*

# clang++ -Wall -march=native -std=c++2a -c chap4.cc -o chap4.o
# clang++ -Wall -march=native -std=c++2a -save-temps chap4.o baseline.cc -o baseline