------------------------------------------------------------------------
LIBERI (LIBrary for numerical evaluation of Electron Repulsion Integral)

Masayuki TOYODA (m-toyoda@jaist.ac.jp), Dec/10/2009
------------------------------------------------------------------------

LIBERI is a C library for numerical evaluation of four-center electron
repulsion integrals (ERIs) based on successive reduction of integral
dimension by using Fourier transform.

For details of the computational method, please refer our article [1] 
as well as the CPC article associated with the package:
  

------------------------------------------------------------------------
FILES IN ARCHIVE:

./README.txt        : this document
./makefile          : sample makefile 
./source/           : source files of LIBERI
./test/             : source files, shell scripts and output examples 
                      of demo programs


------------------------------------------------------------------------
SOURCE FILES:

A brief description of the source files is shown in the list below. 
Those subroutines described in the text and Fig. 3 of our CPC article 
are also shown in the parenthesis.

./source/eri.c 
    High-level functions (ERI_Overlap, ERI_Integral)

./source/eri_ll.c
    Low-level functions (ERI_Init, ERI_Free, ERI_Size_of_*, 
    ERI_Transform_Orbital, ERI_LL_Gamma, ERI_LL_Alpha, ERI_LL_Overlap, 
    ERI_Transform_Overlap, ERI_GL_Interpolate, ERI_Integral_GL)

./source/eri_gaunt_table.c
    Efficient tabulation of the Gaunt coefficients.

./source/eri_interpolate.c
    Numerical interpolation

./source/eri_sf.c
    Calculation of the special functions such as the spherical Bessel 
    functions and the Gaunt coefficients.

./source/sbt/eri_sbt.c
    Interface for the FSBT methods.

./source/sbt/linear/eri_linfsbt.c
    The linear-FSBT method [2].

./source/sbt/log/eri_fsbt.c
    Siegman-Talman FSBT method [3, 4]

./source/sbt/log/eri_logfsbt.c
    The Log-FSBT method.
   

------------------------------------------------------------------------
HOW TO COMPILE:

A sample makefile is supplied with the package of LIBERI. Before 
compiling the archive library, liberi.a, you have to edit the following 
section of the makefile 

  CC   = gcc
  OPT  = -O3 -w -static 

so that CC specifies your C compiler and OPT holds the optimization 
options for the compiler. In addition, you may have to edit the following
section of the makefile
 
  FFTW3_INCLUDE_DIR = $HOME/local/include

where FFTW3_INCLUDE_DIR should point to the directory where fftw3.h is
located. This macro is needed only when the directory of FFTW is not 
specified in other macros such as INCLUDE.

After you edit the makefile, you can compile liberi.a by typing
 
  $ make


------------------------------------------------------------------------
HOW TO LINK LIBERI WITH YOUR PROGRAM:

You can combine LIBERI with your program by including the header file
"eri.h" and by linking it with the archive file "liberi.a". 

Also, you have to link your program with two other external libraries,
namely BLAS and FFTW3. For the installation of them, please access the 
following sites:

BLAS: http://www.betlib.org/blas/
FFTW: http://www.fftw.org/


------------------------------------------------------------------------
DEMO AND SELF-CHECK PROGRAM:

Source files for a simple demo program and for a comprehensive 
self-check program are also supplied in the package. Before you use 
LIBERI in your program, you should perform the self-check test and 
verify that LIBERI can compute ERIs with sufficient accuracy on your 
computers. The source files will be useful for you as a practical 
tutorial how to use the functions of LIBERI.

To compile the programs, you have to edit the following section of the 
makefile 

  LIB  = -lm -g2c -L$(HOME)/local/lib -lfftw3 -lcblas -llapack

so that the external libraries (BLAS and FFTW3) are properly linked 
with the demo programs. In the makefile accompanied with the package, 
several examples of LIB for different computational environments are 
included.

If you have a trouble in linking with the BLAS library, you might want 
to try using FORTRAN compiler to load it. In that case, your FORTRAN 
compiler should be specified by the LOADER variable. Note that FORTRAN 
compiler is not always necessary because LIBERI contains no FORTRAN 
source code. In our experience, LOADER is needed only when we try to 
compile the demo programs with using gcc.

Then, by typing

  $ make demo self-check

or simply
 
  $ make all

a simple demo program ./test/demo/liberi-demo and a comprehensive 
self-check program ./test/selfcheck/liberi-selfcheck will be created.

Shell scripts to run those programs are also provided under ./test 
directory. For example, you can run the demo program by  

  $ cd test
  $ sh run_demo.sh

This program calculates a single ERI for GTO basis functions, which 
will be completed within a second in modern computers.

A typical output of the demo program is as follows:

  $ sh run_demo.sh
    0.007609   -0.000000

where the first value is the real part of the integral and the second
is the imaginary part, both in Hartree units. The result should be 
enough close to the analytical value (0.007608 Hartree).

A comprehensive self-check test can be run by

  $ sh run_selfcheck.sh

This program performs calculations of (1) ERIs for analytical GTO 
basis functions including p and d orbitals, (2) ERIs for numerical PAO 
basis functions, (3) derivatives of ERIs, and (4) ERIs with screening 
of the Coulomb interaction. This will take several minutes.

An example of output is provided in ./test/out_selfcheck. At the bottom
of the output, the results of the calculations are summarized, where the  
numerical errors are evaluated by comparing the result with analytical 
values or numerical reference values. If all the numerical errors are 
less than a threshold (1e-5 Hartree), the following line will be
appeared at the end of the output:
 
  SELF-CHECK TEST HAS COMPLETED SUCCESSFULLY.

If the self-check failed, we kindly ask you to contact us, in order to 
improve the portability of LIBERI, in the following e-mail address: 
 
  Corresponding author: Dr. Masayuki Toyoda, m-toyoda@jaist.ac.jp 

by supplying the output, error messages if any, the version of LIBERI
you used (the digits in the package name), and information of your 
computational environment (operating system, compilers, versions of 
BLAS and FFTW, and so on). We appreciate your cooperation in advance.



------------------------------------------------------------------------
REFERENCES:

[1] M. Toyoda and T. Ozaki, 
    "Numerical evaluation of electron repulsion integrals for 
    pseudoatomic orbitals and their derivatives,"
    The Journal of Chemical Physical, Vol. 130, pp. 124114 (2009).

[2] M. Toyoda and T. Ozaki,
    "Fast spherical Bessel transform via Fast Fourier transform and 
    recurrence formula",
    submitted to Computer Physics Communications.

[3] A.E. Siegman, 
    "Quasi fast Hankel transform,"
    Optics Letters, Vol. 1, pp. 13 (1977).

[4] J.D. Talman,
    "Numerical Fourier and Bessel transform in logarithmic variables,"
    Journal of Computational Physics, Vol. 29, pp. 35 (1978).
