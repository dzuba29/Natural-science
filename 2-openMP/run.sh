#!/bin/bash

omp_threads=2
t=1
n=1000
X=15
Y=15

export OMP_NUM_THREADS=$omp_threads
./main -t $t -n $n -x $X -y $Y -l
