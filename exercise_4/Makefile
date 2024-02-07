CC = clang++
FLAGS = -std=c++2a
OPT_FLAG? = O0

p1: ex4
	${CC} ${FLAGS} -${OPT_FLAG} -o ex4_p1_a_${OPT_FLAG} ex4_p1_a.cc ex4.o
	${CC} ${FLAGS} -${OPT_FLAG} -o ex4_p1_b_${OPT_FLAG} ex4_p1_b.cc ex4.o
p2: ex4
	${CC} ${FLAGS} -${OPT_FLAG} -o ex4_p2_a_${OPT_FLAG} ex4_p2_a.cc ex4.o
	${CC} ${FLAGS} -${OPT_FLAG} -o ex4_p2_b_${OPT_FLAG} ex4_p2_b.cc ex4.o
	${CC} ${FLAGS} -${OPT_FLAG} -o ex4_p2_c_${OPT_FLAG} ex4_p2_c.cc ex4.o

ex4:
	${CC} ${FLAGS} -${OPT_FLAG} -c ex4.cc

clean:
	rm -f ex4.o ex4_p1_a_O0 ex4_p1_a_O3 ex4_p1_b_O0 ex4_p1_b_O3 \
		ex4_p2_a_O0 ex4_p2_a_O3 ex4_p2_b_O0 ex4_p2_b_O3 ex4_p2_c_O0 ex4_p2_c_O3