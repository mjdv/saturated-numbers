# saturated-numbers

A saturated number is a number n such that for all primes p, p divides n if and only if p - 1 divides n.

The single Python file in this repository lists all saturated numbers up to some cutoff, set by default to 1e1000.

The code skips verifying that 619458 does not divide a saturated number below the cutoff, since it takes up the vast majority of the computation time. This was verified (by me) up to a cutoff of 1e1000. If you wish to search up to a larger cutoff, or verify that no such number exists yourself, initialize the variable `exceeded_bases` to the empty list.
