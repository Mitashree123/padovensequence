#!/usr/bin/env python3
"""
primitive_divisors.py

Find all indices m in [4, N] for which the Padovan number P_m has
NO primitive prime divisor.

Padovan sequence used throughout:
    P_0 = 0, P_1 = P_2 = 1,  P_{n+3} = P_{n+1} + P_n.

A prime p is a *primitive divisor* of P_m if p divides P_m but p does
not divide any earlier term P_j (1 <= j < m).

Method (no factorisation is needed):
    Start with d = P_m. For j = 4, ..., m-1, divide d by gcd(d, P_j)
    repeatedly until gcd(d, P_j) = 1. At the end, d is the largest divisor
    of P_m that is coprime to P_1 P_2 ... P_{m-1}. Hence P_m has a
    primitive divisor if and only if d > 1.
    (P_1 = P_2 = P_3 = 1, so starting at j = 4 loses nothing.)

Usage:
    python3 primitive_divisors.py            # N = 1000 (as in the paper)
    python3 primitive_divisors.py 2000       # any other bound

Expected output for N = 1000:
    [5, 7, 10, 11, 12, 13, 14, 16, 21, 23, 32, 45, 58, 61]

Requires only the Python standard library.
"""

import sys
from math import gcd


def padovan(n_max):
    """Return the list [P_0, P_1, ..., P_{n_max}] as exact integers."""
    P = [0, 1, 1]
    for n in range(3, n_max + 1):
        P.append(P[n - 2] + P[n - 3])
    return P


def has_primitive_divisor(m, P):
    """True if P[m] has a prime factor dividing no P[j] with 1 <= j < m."""
    d = P[m]
    for j in range(4, m):
        g = gcd(d, P[j])
        while g > 1:            # remove every prime that P[m] shares with P[j]
            d //= g
            g = gcd(d, P[j])
        if d == 1:              # nothing left: no primitive divisor
            return False
    return d > 1


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    P = padovan(N)
    print(f"P_{N} has {len(str(P[N]))} decimal digits.")
    exceptions = [m for m in range(4, N + 1) if not has_primitive_divisor(m, P)]
    print(f"Indices m in [4, {N}] such that P_m has no primitive divisor:")
    print(exceptions)


if __name__ == "__main__":
    main()
