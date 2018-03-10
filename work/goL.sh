#!/bin/bash --norc
#
#$-S /bin/bash  # interpreting shell
#$-q all.q
#$-cwd        # execute for the current working directory
#$-j y        # y/n: merge stderr into stdout
#$-pe mpi 10   # total num slots. NSLOTS is set to N of PE ; <-- CHECK!!

#OMP_NUM_THREADS=1  # 1 thread per process
#export OMP_NUM_THREADS

mpirun openmxgmo -runtestL -nt 1
