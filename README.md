# Smallest-integer-with-n-distinct-prime-factors
Searches an interval for the smallest integer with at least n distinct prime factors.

This program searches for the smallest integer x with a <= x <= b with x having at least n distinct prime factors.

Author: Simon Goater Dec 2025

Usage: ./smallestwithmprimefactors3.bin numthreads a b n 

1 <= a <= b < 2^64

1 <= n <= 30

E.g. ./smallestwithmprimefactors3.bin 4  10000000000000000000 18000000000000000000 13

10000000000255252260 is the smallest integer in [10000000000000000000, 18000000000000000000] with at least 13 distinct prime factors.

