# Computational methods of linear algebra

This repository contains the solutions of the laboratory work offered as
the part of the "Computational methods of algebra" course at the Faculty of
Applied Mathematics and Computer Science (BSU, 2019).

The base of the solution is the template Matrix class. It can be used for
storing matrices with floating-point values. It includes implementations of
several linear algebra algorithms, in particular:
* Arithmetic operators (including counting the product of matrices).
* TLU decomposition (T stands for the rows permutations) for any square 
non-singular matrix.
* LDL decomposition for square symmetric matriсes.
* Solving the systems of the linear equations using TLU / LDL decomposition.
* Counting the inverse matrix using TLU decomposition.

Some algorithms for the special types of matrices are implemented, for
example, the optimized version of the TLU decomposition algorithm for the
almost triangular matrices.

C++ multithreading features are used inside some algorithms (to efficiently
use all the computer resources). 
[ThreadPool library by Robert Vaser](https://github.com/rvaser/thread_pool) is
sometimes responsible for managing the threads.

More detailed descriprion of the algorithms can be found inside source
files (as comments) and inside the reports, which are available in the form
of Jupyter Notebooks (the Russian language is used there).
