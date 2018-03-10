# -*- coding:utf-8; mode:shell-script; -*-
==== FOR TEST RUN ====

$ cd work  # change directory to the work of openmx

# For any PC
$ mpirun -np 1 openmxgmo Methane.dat -nt 1 > met.std
$ diff ../work.org/input_example/Methane.out met.out |less

# For ivy
$ mpirun -np 1 openmxgmo -runtest -nt 2
$ mpirun -np 2 openmxgmo -runtestL -nt 1  # this job was omitted because of poor CPU

# For hydrangea
$ mpirun -np 2 openmxgmo -runtest -nt 4
$ mpirun -np 8 openmxgmo -runtestL -nt 1

# For spiraea
$ mpirun -np 4 openmxgmo -runtest -nt 1  # use a compute-0-2
$ qsub goL.sh  # use a SGE queueing system

note:
- check the diff output of runtest.result and runtestL.result files
  (any other value of outputs was not checked).
- option -testrun  which probably reads input_example/*.dat
- option -testrunL which probably reads large_example/*.dat
