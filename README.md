# Computational methods of linear algebra

This repository contains the solutions of the laboratory works offered as
the part of the "Computational methods of algebra" course at the Faculty of
Applied Mathematics and Computer Science (BSU, 2019).

The base of the solution is the template Matrix class. It can be used for
storing matrices with floating-point values. It includes implementations of
several linear algebra algorithms, in particular:
* Arithmetic operators (including counting the product of matrices).
* TLU decomposition (T stands for the rows permutations) for any square 
non-singular matrix.
* LDL decomposition for square symmetric matrices.
* Solving the systems of the linear equations using TLU / LDL decomposition.
* Counting the inverse matrix using TLU decomposition.
* Searching for the eigenvalues and eigenvectors with help of several
algorithms (QR algorithm, Power iteration algorithm, Danilevsky algorithm).

Some algorithms for the special types of matrices are implemented, for
example, an optimized version of the TLU decomposition algorithm for the
almost triangular matrices or an algorithm for solving the systems of
linear equations with the tridiagonal matrices.

Template Polynomial class is implemented to work with the characteristic
polynomials of the matrices. It includes arithmetic and differentiation
operators as well as some methods to find the real roots of the polynomial
(Newton and bisection algorithms are used for that).

C++11 multithreading features are used inside some algorithms.
[ThreadPool library by Google](
https://github.com/google/or-tools/blob/v7.4/ortools/base/threadpool.h)
is sometimes responsible for managing the threads.

More detailed descriptions of the algorithms can be found inside source
files (as comments) and inside the reports, which are available in the form
of Jupyter Notebooks (the Russian language is used there).

## Notes

* The project was created and tested on Linux-based systems. There are no
guarantees that it will work correctly in Windows (in particular,
multithreading features may work incorrectly when using the MinGW compiler 
on Windows).

* Some C++2a features are used inside the project (for example, concepts), so
you need some modern compiler to build it. clang++-10 was used to create
this solution (with CLion IDE, which supports C++2a concepts feature), so
it's recommended to work with it (but other compilers can still be used, if
they support necessary features).

