CC=g++
FLAGS = -fopenmp 
#  -O3
parallel: clean parallel.cc
	${CC} parallel.cc -o bin/parallel ${FLAGS}

original: md1d.cpp
	${CC} md1d.cpp -o bin/md1d

clean:
	rm -rf *.dSYM bin/*
