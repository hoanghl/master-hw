FLAGS=-g -O3 -march=native -std=c++17 
COMPILER=clang++


c2v0: c2v0.cc
	${COMPILER} ${FLAGS} c2v0.cc chap2.cc -o c2v0

write: write.cc
	${COMPILER} ${FLAGS} write.cc chap2.cc -o write

chap2.o: chap2.cc
	${COMPILER} ${FLAGS} -c chap2.cc

clean:
	rm -f chap2.o c2v0.o c2v0 write