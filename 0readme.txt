# -*- coding:utf-8; mode:shell-script; -*-
2017/11/17
openmx3.8 + patch3.8.3 をコンパイルし直した。

テストラン
( https://orange.eit.hirosaki-u.ac.jp/nc2/html/htdocs/?page_id=197 )

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

check the diff output of runtest.result and runtestL.result files
(any other value of outputs was not checked).

note:
-testrun  which probably reads input_example/*.dat
-testrunL which probably reads large_example/*.dat
