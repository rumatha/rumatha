"""
Numbers functions.
"""

import sortedcontainers

#===================================================================================================
# Pythagorian triplets.
#===================================================================================================

def pythagorian_triplets(max_n, max_k):
    """
    Sorted set of pythagorian triplets.

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

if __name__ == '__main__':
    print(unique_quadratic_residues(29, True))
    print(pythagorian_triplets(2, 5))

#===================================================================================================