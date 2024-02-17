CC=gfortran

all: p2_a p2_b p2_c p2_d p2_e p2_f

p2_a: omit_frame_pointer.f90
	${CC} -O0 -fno-omit-frame-pointer -o p2_a omit_frame_pointer.f90
p2_b: omit_frame_pointer.f90
	${CC} -O0 -fomit-frame-pointer -o p2_b omit_frame_pointer.f90
p2_c: omit_frame_pointer.f90
	${CC} -O1 -fno-omit-frame-pointer -o p2_c omit_frame_pointer.f90
p2_d: omit_frame_pointer.f90
	${CC} -O1 -fomit-frame-pointer -o p2_d omit_frame_pointer.f90
p2_e: omit_frame_pointer.f90
	${CC} -O2 -fno-omit-frame-pointer -o p2_e omit_frame_pointer.f90
p2_f: omit_frame_pointer.f90
	${CC} -O2 -fomit-frame-pointer -o p2_f omit_frame_pointer.f90
clean:
	rm -f p2_a p2_b p2_c p2_d p2_e p2_f funcmod.mod