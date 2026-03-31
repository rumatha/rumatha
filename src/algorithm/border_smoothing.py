"""
Border smoothing.
Smooth border between two domains of surface unstructured mesh.
"""

import sys
import numpy as np

#===================================================================================================

class Template:
    """
    Template of smoothing.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, lo, hi, ubaln, bbaln):
        """
        Constructor.

        Parameters
        ----------
        lo : int
            Lo point number.
        hi : int
            Hi point number.
        ubaln : int
            Balance change of U domain if template is applied.
        bbaln : int
            Border balance change if template is applied.
        """

        if hi - lo <= 1:
            raise Exception(f'impossible template {lo}, {hi}')
        if bbaln != -1:
            raise Exception('template with border balance != -1 is not supported')

        self.lo = lo
        self.hi = hi
        self.ubaln = ubaln
        self.dbaln = -self.ubaln
        self.bbaln = bbaln

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'T[{self.lo}, {self.hi}](U:{self.ubaln},B:{self.bbaln})'

    #-----------------------------------------------------------------------------------------------

    def has_conflict_with(self, t):
        """
        Check if template has conflict with another template.

        Parameters
        ----------
        t : Template
            Another template.

        Returns
        -------
        bool
            True - if there is conflict,
            False - otherwise.
        """

        if (self.lo >= t.hi) or (t.lo >= self.hi):
            return False

        return True

#===================================================================================================

class Templates:
    """
    Templates class.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, ts=[]):
        """
        Constructor.
        """

        self.items = []
        self.add(ts)

    #-----------------------------------------------------------------------------------------------

    def add(self, t):
        """
        Add template of templates.

        Parameters
        ----------
        t : Template | [Template]
            Template or list of templates to add.
        """

        if isinstance(t, Template):
            self.items.append(t)
        elif isinstance(t, list):
            for ti in t:
                self.add(ti)
        else:
            raise Exception('can not add object to templates')

    #-----------------------------------------------------------------------------------------------

    def count(self):
        """
        Count of elements.

        Returns
        -------
        int
            Count of elements.
        """

        return len(self.items)

    #-----------------------------------------------------------------------------------------------

    def min_ubaln(self):
        """
        Minimum possible U balance.

        Returns
        -------
        int
            Minimum possible U balance.
        """

        return sum([t.ubaln for t in self.items if t.ubaln < 0])

    #-----------------------------------------------------------------------------------------------

    def max_ubaln(self):
        """
        Maximum possible U balance.

        Returns
        -------
        int
            Maximum possible U balance.
        """

        return sum([t.ubaln for t in self.items if t.ubaln > 0])

#===================================================================================================

def build_b_matrix(ts, is_log=False):
    """
    Build B matrix for templates.

    Parameters
    ----------
    ts : Templates
        Templates.
    is_log : bool
        Log actions.

    Returns
    -------
    np.array
        3-dimensional array.
    """

    # Allocate memory.
    # First dimension - count of templates.
    # Second dimension - possible variation of balance between domains.
    # Third dimension - 2 (use or do not use current template).
    tn = ts.count()
    min_ubaln, max_ubaln = ts.min_ubaln(), ts.max_ubaln()
    b = np.full((tn, max_ubaln - min_ubaln + 1, 2), sys.maxsize, 'int')

    # Init for last template.
    ci = tn - 1
    ct = ts.items[ci]
    b[ci][0][0] = 0
    b[ci][ct.ubaln][1] = -1

    # Process all templates from last to first.
    while True:

        # Shift templates.
        ni, nt = ci, ct
        ci = ci - 1
        if ci < 0:
            break
        ct = ts.items[ci]

        # ct - current template
        # nt - next template (already processed)

        # I. Ignore current template.
        for j in range(min_ubaln, max_ubaln + 1):
            b[ci][j][0] = min(b[ni][j][0], b[ni][j][1])

        if not ct.has_conflict_with(nt):

            # II. Check for no conflict.
            for j in range(min_ubaln, max_ubaln + 1):
                if b[ci][j][0] != sys.maxsize:
                    b[ci][j + ct.ubaln][1] = b[ci][j][0] - 1
        else:

            # III. Check for conflict.
            for j in range(min_ubaln, max_ubaln + 1):
                if b[ni][j][0] != sys.maxsize:
                    b[ci][j + ct.ubaln][1] = b[ni][j][0] - 1

    # Print matrix b.
    if is_log:
        def fstr(v):
            return v if v != sys.maxsize else 'F'
        for i in range(tn):
            str = ' '.join([f'[{fstr(b[i][j][0])} {fstr(b[i][j][1])}]'
                            for j in range(min_ubaln, max_ubaln + 1)])
            print(str)

    return b

#===================================================================================================

if __name__ == '__main__':

    # Test border smoothing.
    ts = Templates([Template( 0,  2, 1, -1),
                    Template( 1,  4, -3, -1),
                    Template( 3,  5, 1, -1),
                    Template( 5,  7, -1, -1),
                    Template( 6,  8, 1, -1),
                    Template( 7, 10, -3, -1),
                    Template( 9, 11, 1, -1),
                    Template(10, 12, -1, -1),
                    Template(11, 14, +3, -1),
                    Template(13, 15, -1, -1),
                    Template(14, 16, +1, -1)])
    build_b_matrix(ts, True)

#===================================================================================================
