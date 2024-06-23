# CC = g++
# CC_GPU = nvcc

# BIN = bin
# BUILD = build
# SRC = src
# DIR_LECTURE := lecture
# DIR_HW = hw

# FLAGS = -march=native -std=c++2a  -fopenmp -mavx -O3

# hw: clean chap4.o ppc_kernels.o ${DIR_HW}/main_gpu.cu
# 	${CC_GPU} chap4.o ppc_kernels.o ${DIR_HW}/main_gpu.cu -o ${BIN}/run_hw_gpu

# ppc_kernels.o: ${DIR_HW}/ppc_kernels.cu
# 	${CC_GPU} -c ${DIR_HW}/ppc_kernels.cu

# lec: clean chap4.o ${DIR_LECTURE}/gpu_v0.cu
# 	${CC_GPU} chap4.o ${DIR_LECTURE}/gpu_v0.cu -o bin/gpu

# lecture: clean chap4.o ${DIR_LECTURE}/lec_cpu.cc ${DIR_LECTURE}/lec_gpu.cu ${DIR_LECTURE}/lec_gpu_v3.cu
# 	${CC} chap4.o ${DIR_LECTURE}/lec_cpu.cc -o ${BIN}/run_lec_cpu -g ${FLAGS}
# 	${CC_GPU} chap4.o ${DIR_LECTURE}/lec_gpu.cu -o ${BIN}/run_lec_gpu
# 	${CC_GPU} chap4.o ${DIR_LECTURE}/lec_gpu_v3.cu -o ${BIN}/run_lec_gpu_3

# chap4: 
# 	${CC} ${FLAGS} -I ./include/ -c ${SRC}/cpu.cc -o ${BUILD}/cpu.o
# 	${CC} ${FLAGS} -I ./include/ -c ${SRC}/chap4.cc -o ${BUILD}/chap4.o

# clean:
# 	rm -f ${BUILD}/* ${BIN}/*

CC = g++
CC_GPU = nvcc

BIN = bin
BUILD = build
SRC = src
FLAGS = -march=native -std=c++2a -fopenmp -mavx -O3

main: utils gpu cpu
	${CC_GPU} -I ./include/ ${BUILD}/gpu.o ${BUILD}/utils.o ${SRC}/main_gpu.cu -o ${BIN}/gpu
	${CC} ${FLAGS} -I ./include/ ${BUILD}/cpu.o ${BUILD}/utils.o ${SRC}/main_cpu.cc -o ${BIN}/cpu

gpu: 
	${CC_GPU} -I ./include/ -c ${SRC}/gpu.cu -o ${BUILD}/gpu.o
cpu: 
	${CC} ${FLAGS} -I ./include/ -c ${SRC}/cpu.cc -o ${BUILD}/cpu.o

utils:
	${CC} ${FLAGS} -I ./include/ -c ${SRC}/utils.cc -o ${BUILD}/utils.o

clean:
	rm -f ${BUILD}/* ${BIN}/*