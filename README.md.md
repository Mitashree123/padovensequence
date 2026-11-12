# Multiplicatively dependent tuples of Padovan numbers — data and code

This repository contains the data and the programs for the paper

> [P. K. BHOI, M. BEHERA, AND M. BEHERA], *Multiplicatively dependent tuples of Padovan numbers*.

## The problem

The Padovan sequence is

    P_0 = 0,  P_1 = P_2 = 1,  P_{n+3} = P_{n+1} + P_n,

that is, 0, 1, 1, 1, 2, 2, 3, 4, 5, 7, 9, 12, 16, 21, 28, … (OEIS A000931 shifted by three places).

A tuple (P_{m_1}, …, P_{m_k}) with 4 ≤ m_1 < … < m_k is **qualifying** if the numbers are multiplicatively dependent,

    P_{m_1}^{y_1} · P_{m_2}^{y_2} ⋯ P_{m_k}^{y_k} = 1   (y_i integers, not all zero),

while any k − 1 of them are multiplicatively independent. The vector y = (y_1, …, y_k) is then unique up to sign; we take it primitive with y_1 > 0. Put

    S = y_1 + … + y_k,    T = m_1 y_1 + … + m_k y_k.

The tuple is **degenerate** if S = T = 0 and **non-degenerate** otherwise. (In the paper this is shown to be equivalent to whether ∏ a^{-y_i} α^{-m_i y_i} equals 1, where P_n = aα^n + bβ^n + cγ^n.)

## Contents

| File | Description |
|---|---|
| `padovan_tuples.csv` | All 1993 qualifying tuples with m_k ≤ 1000 |
| `primitive_divisors.py` | Finds the indices m ≤ N for which P_m has no primitive prime divisor |
| `enumerate_tuples.py` | Enumerates all qualifying tuples and writes `padovan_tuples.csv` |
| `requirements.txt` | Python dependencies |

## The data file `padovan_tuples.csv`

One row per qualifying tuple, sorted by k and then by the indices.

| Column | Meaning |
|---|---|
| `k` | number of terms in the tuple |
| `indices` | m_1 < m_2 < … < m_k, separated by spaces |
| `exponent_vector_y` | the relation vector y_1 … y_k, separated by spaces |
| `S=sum_y` | S = Σ y_i |
| `T=sum_m*y` | T = Σ m_i y_i |
| `class` | `degenerate` (S = T = 0) or `non-degenerate` |
| `identity` | the relation written with positive exponents on both sides, e.g. `P10 * P12 = P11^2` |
| `common_value` | the integer that both sides of the identity equal |

Example rows:

    k,indices,exponent_vector_y,S=sum_y,T=sum_m*y,class,identity,common_value
    2,4 5,1 -1,0,-1,non-degenerate,P4 = P5,2
    3,10 11 12,1 -2 1,0,0,degenerate,P10 * P12 = P11^2,144

Summary:

| k | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | total |
|---|---|---|---|---|---|---|---|---|---|
| degenerate | 0 | 3 | 9 | 1 | 20 | 0 | 1 | 5 | **39** |
| non-degenerate | 8 | 21 | 35 | 189 | 687 | 4 | 84 | 926 | **1954** |
| all | 8 | 24 | 44 | 190 | 707 | 4 | 85 | 931 | **1993** |

## How the list is obtained

1. **Primitive divisors.** In a qualifying tuple, the largest term P_{m_k} has no primitive prime divisor (a prime dividing P_{m_k} but no earlier P_j); otherwise that prime could not cancel in the relation. `primitive_divisors.py` shows, using only gcd computations, that for 4 ≤ m ≤ 1000 this happens only for

       m = 5, 7, 10, 11, 12, 13, 14, 16, 21, 23, 32, 45, 58, 61.

   Hence every qualifying tuple with m_k ≤ 1000 has all its indices in [4, 61].

2. **Enumeration.** `enumerate_tuples.py` factors P_4, …, P_61 (all below 1.6·10^7). A set of indices is a qualifying tuple exactly when the prime-exponent vectors of the corresponding P_m form a *circuit* (they are linearly dependent, but every proper subset is independent). The program finds all circuits, computes the relation vector y exactly, checks each relation in exact integer arithmetic, and classifies the tuple using S and T.

**Scope.** The list is complete for m_k ≤ 1000. It is the complete list of all qualifying tuples, for every k, if P_m has a primitive divisor for every m ≥ 62 (Conjecture 9.2 of the paper).

## Running the code

Requirements: Python 3.8 or later and `sympy` (`pip install -r requirements.txt`).

    python3 primitive_divisors.py          # takes a few seconds
    python3 enumerate_tuples.py            # takes a few minutes; writes padovan_tuples.csv

`primitive_divisors.py` accepts an optional bound, e.g. `python3 primitive_divisors.py 2000`.

## How to cite

If you use these data or programs, please cite the paper and this archive:

> [Author name], *Multiplicatively dependent tuples of Padovan numbers: data and code*, version 1.0, Zenodo (2026), doi:10.5281/zenodo.XXXXXXX.

## Licence

The data (`padovan_tuples.csv`) are released under the Creative Commons Attribution 4.0 International licence (CC BY 4.0), see `LICENSE-DATA`. The code is released under the MIT licence, see `LICENSE`.
