# Instructions

## 1. Prerequisites

All codes for each problem share the same header `ex4.hpp`.
Please replace variable `CC` in `Makefile` to `g++`.

## 2. Make

`Make` command has following format:
`make OPT_FLAG=<type of optimization flag: O0 or O3> <target>`

Example:

`make OPT_FLAG=O3 p1`

## 3. For problem 3

First, run the code with command:
`    python ex4_p3/ex4_p3_gencode.py`
It will generate **csv** file containing the running time. Then use file `visualize.ipynb` to plot the figure.
