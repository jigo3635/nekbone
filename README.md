# Nekbone with LIBXSMM/SIMD

All the compilers/flags and variables can be changed in test/example[1-3]/makenek

### General information

Set the F77/CC compilers

```
# Fortran compiler
F77="mpiifort"

# C compiler
CC="mpiicc"
```

### Using the in-hourse SIMD implementation

uncomment out the variable "IFSIMD"

```
# Enable SIMD (default false)
#IFSIMD="true"
```

### Using the libxsmm

uncomment out the variable "IFXSMM"

```
# Enable XSMM (default false)
#IFXSMM="true"
```
In the case, include paths and linked libraries can be set

```
USR_INCDIR="..."
USR_LDFLAGS="..."
```

### file test/SIZE
the polynomial order and the number of elements per core can be changed in the file SIZE

### Compiling and running the code

```
$ cd test
$ ./makenek clean
$ ./makenek test
$ mpirun -n  ./nekbone
```




