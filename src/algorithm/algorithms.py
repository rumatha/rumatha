"""
Algorithms for various purposes.
"""

import sortedcontainers

#===================================================================================================
# Get expressions set.
#===================================================================================================

def get_expressions_set_from_digits(ds):
    """
    Get all expressions from digits.
    Possible operations: (), +, -, *, /, ^.
    We do not leave Z numbers.
    We do not use too large numbers (1^n is available, m^n for n < power_limit is available).

    Parameters
    ----------
    ds : str
        Digits written in string.

    Returns
    -------
    set()
        Set of expressions - tuples (str, str).
    """

    # Limit for power operation.
    power_limit = 10

    n, s = len(ds), sortedcontainers.SortedSet([(int(ds), f'{ds}')])
    if n > 1:
        for i in range(1, n):
            los = get_expressions_set_from_digits(ds[:i])
            his = get_expressions_set_from_digits(ds[i:])
            for lov, loe in los:
                for hiv, hie in his:

                    # Add operation.
                    s.add((lov + hiv, f'({loe} + {hie})'))

                    # Sub operation.
                    s.add((lov - hiv, f'({loe} - {hie})'))

                    # Mul operation.
                    s.add((lov * hiv, f'({loe} * {hie})'))

                    # Div operation (we do not leave Z numbers).
                    if (hiv != 0) and (lov % hiv == 0):
                        s.add((lov // hiv, f'({loe} / {hie})'))

                    # Power operation (we do not leave Z numbers and do not use too large numbers).
                    if hiv > 0:
                        if (lov == 1) or (hiv < power_limit):
                            s.add((lov**hiv, f'({loe}^{hie})'))
    return s

#===================================================================================================

if __name__ == '__main__':
    pass

#===================================================================================================
