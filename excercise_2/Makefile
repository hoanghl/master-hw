CC=gcc

p1:
	${CC} -o problem1 problem_1.c 
p2:
	${CC} -g -o problem2 problem_1.c 
p3:
	${CC} -lm -o problem3 problem_3.c

p4: 
	${CC} -O0 -pg -g -static -o problem4 simpleprofiling/simpleprofiling.c -lm

clean:
	rm -f problem1 problem2 problem3 problem4
