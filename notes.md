# Notes

# 1. CPU parallelism

## Technique 1: Compiling optimization

Normally, compiling Cpp program is done with g++ command. However, there are flags we can use to boost the efficiency.

| Flag           | Meaning                                                                                                                                                                         |
| -------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `march=native` | Compile program to native CPU architecture. This implies the compiled program is tailored for the current CPU architecture.                                                     |
| `-O3`          | A optimizatio flag. This will optimize deeply the program but the downside is that the compiler will add more code to the final compiled program. Furthermore, this is buggier. |

## Technique 2: Linear reading

Modern CPU c is excel as predicting which memory cells are accessed in the future and it will load values in consecutive cells to register beforehand. If our code can take advantage of this predicting mechanism, we will have another method for optimization, namely _Linear reading_.
Looking at the figure below and imagining the case as flatting the matrix to 1D array, says position `[i, j]` in array is converted to position `[n*i + j]` in 1D array. Given this conversion, it's easy to see that consecutive looping in `i` direction is conserved while that is not true in `j` direction. Consequently, to ensure the consecutiveness in `j` direction and therefore utilize CPU memory overhead prediction, a modification to original memory access is needed.
One of possible modification is _transposing_. In particular, if the given matrix is transposed, looping via `j` direction is now consecutive.
![Array](doc/1.jpg)

## Technique 3: Instruction-level parallelism
