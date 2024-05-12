# 2D Molecular Dynamics Simulation

Huy Hoang Le

## 1. Directory structure

This directory contains the following:

- Implementation of serial and parallel version of `md2d`: file `md2d.cc` for serial version, file `md2d_parallel.cc` for parallel version
- Makefile
- benchmarking code and result

## 2. Run

### 2.1. Build

For serial version, I use `clang++` for building. Please modify to `g++` if your system doen't has `clang`.

Run 2 following commands to build serial and parallel version:

```
make md2d
make md2d_parallel
```

Binary files are stored in directory `bin`.

### 2.2. (Optional) Benchmark

To run the benchmarking code in **Python**, run the notebook `evaluation.ipynb` in folder `benchmarking`
