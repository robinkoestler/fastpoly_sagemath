# A fast class for polynomials in SageMath

## Motivation

FHE needs efficient operations in its native ring $\mathbb Z[X] / \langle X^N+1 \rangle $.
Unfortunately, the reduction mod $X^N+1$ and the automorphism evaluations $X \longmapsto X^5$ for example, are slow if implemented in plain SageMath.
This motivates to lift the NTL class, which is sage-inbuilt, for FHE purposes.

## Explanation

The entire class "Poly" works with objects either from the class sage.libs.ntl.ntl_ZZ_pX or sage.libs.ntl.ntl_ZZX for computations mod p resp. in $\mathbb{Z}$.
It tries to not leave these inbuilt Cython classes, however, sometimes this is impossible to avert.

An element of "Poly" has the attributes:

- modulus: gives the current modulus
- N: is the ring degree as above

The following methods used in FHE are included:

- Fast arithmetic with the operators $+, -, *, **, +=, -=, *=, /, <<, >>$, the last two shifting the polynomials coefficients.
- The modular reduction and modular switching in between $\mathbb{Z}_p[X]$ for $p \in \{0,2,3,4,\dots\}$.
- Fast evaluation of automorphisms:
  - as $X \longmapsto X^5$ by $\mathrm{auto5()}$,
  - as $X \longmapsto X^{-1}$ by $\mathrm{auto\_inverse()}$ and
  - as $X \longmapsto X^{5^j}$ by $\mathrm{auto(j)}$.
- Printing an element always outputs a centered representative with coefficients in the interval $[-Q/2, Q/2)$, which contains $\mathrm{centered\_list()}$.
- Generation of random elements of $R_Q$.
- Norm of a polynomial (max of absolute values)
- Rescaling the polynomial (and modulus) $Q$ by a divisor $q$ of $Q$.
- The fast multiplication also makes use of the $\mod X^N + 1$ reduction in $\mathrm{mod\_quo()}$.
- Several other type checks as coefficient manipulation/extraction, comparison functions, conversion to lists, copying, etc.

## Slow parts

Subroutines, which can be improved, include:

1. Rescaling: Because in SageMath's ntl.ntl_ZZ_pX and ntl.ntl_ZZ_p there is no proper division algorithm, we have to use the _right_shift() method. This only works correctly as a division, if the modulus is a multiple of the divisor. In this case it also is reasonable fast, however it would be of great interest to have this division available for general moduli (without using the slow conversion methods).
2. Switching from $\mathbb{Z}[X]$ to $\mathbb{Z}_p[X]$. This is the only slow (100x slower) modulo conversion, which seems unavoidable.

## Benchmarks

We test our code with the FHE standard $N=2^{15}$ and a modulus $Q = 2^{1000}$ in the file "Test.ipynb".
