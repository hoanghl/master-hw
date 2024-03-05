FLAGS=-O3 -march=native -std=c++2a -fopenmp -mavx

COMPILER=g++

c2v7: c2v7.cc chap2.o
	${COMPILER} -o c2v7 c2v7.cc chap2.o ${FLAGS}

c2v5: c2v5.cc chap2.o
	${COMPILER} -o c2v5 c2v5.cc chap2.o ${FLAGS}

c2v4: c2v4.cc chap2.o
	${COMPILER} -o c2v4 c2v4.cc chap2.o ${FLAGS}

c2v3: c2v3.cc chap2.o
	${COMPILER} -o c2v3 c2v3.cc chap2.o ${FLAGS}

c2v2: c2v2.cc chap2.o
	${COMPILER} c2v2.cc chap2.o -o c2v2 ${FLAGS}

c2v1: c2v1.cc chap2.o
	${COMPILER} c2v1.cc chap2.o -o c2v1 ${FLAGS}
	
c2v0_1: c2v0_1.cc chap2.o
	${COMPILER} -fopenmp c2v0_1.cc chap2.o -o c2v0_1 ${FLAGS}
	
c2v0: c2v0.cc
	${COMPILER} c2v0.cc chap2.cc -o c2v0 ${FLAGS}

write: write.cc
	${COMPILER} write.cc chap2.cc -o write ${FLAGS}

chap2.o: chap2.cc
	${COMPILER} -c chap2.cc ${FLAGS}

clean:
	rm -rf chap2.o c2v0.o c2v0 c2v0_1 c2v1_1 write c2v1 c2v2 c2v2_1 c2v3 c2v4 c2v5 c2v7