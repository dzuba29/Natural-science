#!/bin/bash

omp_threads=4
t=1
n=10000
X=50
Y=50

export OMP_NUM_THREADS=$omp_threads
./main -t $t -n $n -x $X -y $Y -l
