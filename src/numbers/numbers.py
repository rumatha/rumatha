"""
Numbers functions.
"""

import sortedcontainers

#===================================================================================================
# Calculations by module.
#===================================================================================================

def quadratic_residues(p, is_0_included=True):
    """
    All quadratic residues.

    Parameters
    ----------
    p : int
        Module.
    is_0_included : bool
        Flag for including 0 into list.

    Returns
    -------
    [int]
        List of quadratic residues.
    """

    s = 0 if is_0_included else 1

    return [(a * a) % p for a in range(s, p)]

#---------------------------------------------------------------------------------------------------

def unique_quadratic_residues(p, is_0_included=True):
    """
    Unique quadratic residues.

    Parameters
    ----------
    p : int
        Module.
    is_0_included : bool
        Flag for including 0 into set.

    Returns
    -------
    sortedcontainers.SortedSet
        Sorted set of quadratic residues.
    """

    return sortedcontainers.SortedSet(quadratic_residues(p, is_0_included))

#===================================================================================================
# Pythagorean triplets.
#===================================================================================================

def pythagorean_triplets(max_n, max_k):
    """
    Sorted set of pythagorean triplets.

    a = k * (n^2 - m^2)
    b = k * 2 * n * m
    c = k * (n^2 + m^2)

    Parameters
    ----------
    max_n : int
        Maximum value of n.
    max_k : int
        Maximum value of k.

    Returns
    -------
    sortedcontainers.SortedSet
        Set of triplets.
    """

    s = sortedcontainers.SortedSet()

    for k in range(1, max_k + 1):
        for m in range(1, max_n + 1):
            for n in range(m + 1, max_n + 1):
                a = k * (n**2 - m**2)
                b = 2 * k * n * m
                c = k * (n**2 + m**2)
                s.add((a, b, c))

    return s

#---------------------------------------------------------------------------------------------------

def pythagorean_triplets_by_module(p):
    """
    Get all pythagorean triplets by module of p.

    Parameters
    ----------
    p : int
        Module.

    Returns
    -------
    sortedcontainers.SortedSet
        Set of triplets.
    """

    s = sortedcontainers.SortedSet()
    r = range(p)

    for v in r:
        for u in range(v, p):
            for k in r:
                a = (k * (u**2 - v**2)) % p
                b = (k * 2 * u * v) % p
                c = (k * (u**2 + v**2)) % p
                s.add((a, b, c))

    return s

#===================================================================================================

if __name__ == '__main__':

    # Calculations by module.
    print('quadratic_residues(17, False) :', quadratic_residues(17, False))
    print('unique_quadratic_residues(17, True) :',
          unique_quadratic_residues(17, True))

    # Pythagorean triplets.
    print('pythagorean_triplets(2, 5) :', pythagorean_triplets(2, 5))
    print('pythagorean_triplets_by_module(7) :', pythagorean_triplets_by_module(7))

#===================================================================================================
