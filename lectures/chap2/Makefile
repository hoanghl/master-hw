FLAGS=-O3 -march=native -std=c++2a -fopenmp 

COMPILER=g++

c2v5: c2v5.cc chap2.o
	${COMPILER} ${FLAGS} -o c2v5 c2v5.cc chap2.o

c2v4: c2v4.cc chap2.o
	${COMPILER} ${FLAGS} -o c2v4 c2v4.cc chap2.o

c2v3: c2v3.cc chap2.o
	${COMPILER} ${FLAGS} -o c2v3 c2v3.cc chap2.o

c2v2: c2v2.cc chap2.o
	${COMPILER} ${FLAGS} c2v2.cc chap2.o -o c2v2

c2v1: c2v1.cc chap2.o
	${COMPILER} ${FLAGS} c2v1.cc chap2.o -o c2v1
	
c2v0_1: c2v0_1.cc chap2.o
	${COMPILER} ${FLAGS} -fopenmp c2v0_1.cc chap2.o -o c2v0_1
	
c2v0: c2v0.cc
	${COMPILER} ${FLAGS} c2v0.cc chap2.cc -o c2v0

write: write.cc
	${COMPILER} ${FLAGS} write.cc chap2.cc -o write

chap2.o: chap2.cc
	${COMPILER} ${FLAGS} -c chap2.cc

clean:
	rm -rf chap2.o c2v0.o c2v0 c2v0_1 c2v1_1 write c2v1 c2v2 c2v2_1 c2v3 c2v4 c2v5