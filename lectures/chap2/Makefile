FLAGS=-g -O3 -march=native -std=c++2a
COMPILER=g++

c2v2_1: c2v2_1.cc chap2.o
	${COMPILER} ${FLAGS} -fopenmp -c c2v2_1.cc
	${COMPILER} ${FLAGS} -fopenmp c2v2_1.o chap2.o -o c2v2_1

c2v2: c2v2.cc chap2.o
	${COMPILER} ${FLAGS} -c c2v2.cc
	${COMPILER} ${FLAGS} c2v2.o chap2.o -o c2v2

c2v1_1: c2v1_1.cc chap2.o
	${COMPILER} ${FLAGS} -fopenmp -c c2v1_1.cc
	${COMPILER} ${FLAGS} -fopenmp c2v1_1.o chap2.o -o c2v1_1

c2v1: c2v1.cc chap2.o
	${COMPILER} ${FLAGS} -c c2v1.cc
	${COMPILER} ${FLAGS} c2v1.o chap2.o -o c2v1
	
c2v0_1: c2v0_1.cc chap2.o
	${COMPILER} ${FLAGS} -fopenmp -c c2v0_1.cc
	${COMPILER} ${FLAGS} -fopenmp c2v0_1.o chap2.o -o c2v0_1
	
c2v0: c2v0.cc
	${COMPILER} ${FLAGS} c2v0.cc chap2.cc -o c2v0

write: write.cc
	${COMPILER} ${FLAGS} write.cc chap2.cc -o write

chap2.o: chap2.cc
	${COMPILER} ${FLAGS} -c chap2.cc

clean:
	rm -rf chap2.o c2v0.o c2v0 c2v0_1.o c2v0_1 c2v1_1.o c2v1_1 write *.dSYM c2v1.o c2v1 c2v2.o c2v2