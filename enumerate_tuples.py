#!/usr/bin/env python3
"""
enumerate_tuples.py

Enumerate all qualifying tuples of Padovan numbers with largest index <= 61
and write them to padovan_tuples.csv.

Definitions (see the paper):
  * A tuple (P_{m_1}, ..., P_{m_k}), 4 <= m_1 < ... < m_k, is *qualifying*
    if P_{m_1}, ..., P_{m_k} are multiplicatively dependent but any k-1 of
    them are multiplicatively independent.
  * Its relation vector y = (y_1, ..., y_k) is the primitive integer vector
    with  P_{m_1}^{y_1} ... P_{m_k}^{y_k} = 1  (unique up to sign; we take
    y_1 > 0).
  * S = sum(y_i),  T = sum(m_i * y_i).
    The tuple is *degenerate* if S = T = 0, and *non-degenerate* otherwise.

Why only indices <= 61:
  In a qualifying tuple the largest term P_{m_k} has no primitive divisor
  (Lemma 9.1 of the paper). By primitive_divisors.py, the indices m <= 1000
  with this property are all <= 61. Hence every qualifying tuple with
  m_k <= 1000 has all its indices in [4, 61].

Method:
  1. Factor P_4, ..., P_61 (all are below 1.6 * 10^7) and form the vector
     v_m of prime exponents of P_m.
  2. A set of indices is a qualifying tuple exactly when it is a *circuit*
     of the vectors v_m: the vectors are linearly dependent, but every
     proper subset is independent.
  3. Compute the lattice of all relations among v_4, ..., v_61 and keep only
     the indices that occur in some relation (21 indices). Their vectors span
     a space of dimension 8, so a circuit has at most 9 elements.
  4. Test every subset of these 21 indices with at most 10 elements.
     Ranks are computed modulo the prime 2^61 - 1 for speed; the relation
     vector of each circuit is then computed in exact rational arithmetic,
     and the relation is checked in exact integer arithmetic.

Expected result: 1993 tuples (39 degenerate, 1954 non-degenerate).
Running time: a few minutes.

Requires: Python 3 and sympy.
"""

import csv
import itertools
from collections import Counter

import sympy as sp

MAX_INDEX = 61             # largest possible m_k (see docstring)
MAX_SIZE = 10              # test subsets up to this size (rank bound gives 9)
PRIME = (1 << 61) - 1      # modulus for fast rank computations
OUTPUT = "padovan_tuples.csv"


def padovan(n_max):
    """Return [P_0, ..., P_{n_max}]."""
    P = [0, 1, 1]
    for n in range(3, n_max + 1):
        P.append(P[n - 2] + P[n - 3])
    return P


def rank_mod_p(columns):
    """Rank of the matrix with the given integer columns, modulo PRIME."""
    if not columns:
        return 0
    n_rows = len(columns[0])
    A = [[col[r] % PRIME for col in columns] for r in range(n_rows)]
    n_cols = len(columns)
    rank = 0
    for c in range(n_cols):
        pivot = next((i for i in range(rank, n_rows) if A[i][c]), None)
        if pivot is None:
            continue
        A[rank], A[pivot] = A[pivot], A[rank]
        inv = pow(A[rank][c], PRIME - 2, PRIME)
        A[rank] = [x * inv % PRIME for x in A[rank]]
        for i in range(n_rows):
            if i != rank and A[i][c]:
                f = A[i][c]
                A[i] = [(A[i][t] - f * A[rank][t]) % PRIME for t in range(n_cols)]
        rank += 1
    return rank


def primitive_relation(columns):
    """Primitive integer vector y with sum y_i * columns[i] = 0 and y_1 > 0
    (the kernel must be one-dimensional)."""
    A = sp.Matrix(columns).T                 # rows = primes, columns = indices
    kernel = A.nullspace()
    assert len(kernel) == 1, "kernel is not one-dimensional"
    v = kernel[0]
    denom = sp.ilcm(*[sp.Rational(x).q for x in v], 1)
    y = [int(sp.Rational(x) * denom) for x in v]
    g = sp.igcd(*y)
    y = [x // g for x in y]
    if y[0] < 0:
        y = [-x for x in y]
    return y


def format_identity(indices, y):
    """Write P_{m_1}^{y_1} ... = 1 as  (positive part) = (negative part)."""
    def term(m, e):
        return f"P{m}" + (f"^{e}" if e > 1 else "")
    left = [term(m, e) for m, e in zip(indices, y) if e > 0]
    right = [term(m, -e) for m, e in zip(indices, y) if e < 0]
    return " * ".join(left) + " = " + " * ".join(right)


def main():
    P = padovan(MAX_INDEX)
    indices = list(range(4, MAX_INDEX + 1))

    # Step 1: prime-exponent vectors v_m.
    factorisation = {m: sp.factorint(P[m]) for m in indices}
    primes = sorted({p for m in indices for p in factorisation[m]})
    v = {m: [factorisation[m].get(p, 0) for p in primes] for m in indices}

    # Step 3: indices that occur in at least one relation.
    relation_basis = sp.Matrix([v[m] for m in indices]).T.nullspace()
    involved = sorted({indices[i] for b in relation_basis
                       for i in range(len(indices)) if b[i] != 0})
    print(f"Rank of the relation lattice: {len(relation_basis)}")
    print(f"Indices occurring in some relation ({len(involved)}): {involved}")
    print(f"Dimension spanned by their vectors: "
          f"{rank_mod_p([v[m] for m in involved])}")

    # Step 4: all circuits among the involved indices.
    rows = []
    for size in range(2, MAX_SIZE + 1):
        for subset in itertools.combinations(involved, size):
            cols = [v[m] for m in subset]
            if rank_mod_p(cols) == size:
                continue                        # independent: not a relation
            if any(rank_mod_p(cols[:j] + cols[j + 1:]) < size - 1
                   for j in range(size)):
                continue                        # a proper subset is dependent
            y = primitive_relation(cols)

            # Exact check of P_{m_1}^{y_1} ... P_{m_k}^{y_k} = 1.
            left = right = 1
            for m, e in zip(subset, y):
                if e > 0:
                    left *= P[m] ** e
                else:
                    right *= P[m] ** (-e)
            assert left == right, f"relation fails for {subset}"

            S = sum(y)
            T = sum(m * e for m, e in zip(subset, y))
            rows.append([
                size,
                " ".join(map(str, subset)),
                " ".join(map(str, y)),
                S,
                T,
                "degenerate" if S == 0 and T == 0 else "non-degenerate",
                format_identity(subset, y),
                left,
            ])

    # Output.
    with open(OUTPUT, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["k", "indices", "exponent_vector_y", "S=sum_y",
                         "T=sum_m*y", "class", "identity", "common_value"])
        writer.writerows(rows)

    counts = Counter((r[0], r[5]) for r in rows)
    print(f"\nTotal qualifying tuples: {len(rows)}")
    print(f"{'k':>3} {'degenerate':>11} {'non-degenerate':>15}")
    for k in sorted({r[0] for r in rows}):
        print(f"{k:>3} {counts[(k, 'degenerate')]:>11} "
              f"{counts[(k, 'non-degenerate')]:>15}")
    print(f"Written to {OUTPUT}")


if __name__ == "__main__":
    main()
