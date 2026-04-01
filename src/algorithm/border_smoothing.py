"""
Border smoothing.
Smooth border between two domains of surface unstructured mesh.
"""

import sys
import random
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

    def __eq__(self, t):
        """
        Check templates equal.

        Parameters
        ----------
        t : Template
            Another template.

        Returns
        -------
        bool
            True - if templates are equal,
            False - otherwise.
        """

        return (self.lo, self.hi) == (t.lo, t.hi)

    #-----------------------------------------------------------------------------------------------

    def __ne__(self, t):
        """
        Check templates not equal.

        Parameters
        ----------
        t : Template
            Another template.

        Returns
        -------
        bool
            True - if templates are not equal,
            False - otherwise.
        """

        return not (self == t)

    #-----------------------------------------------------------------------------------------------

    def __gt__(self, t):
        """
        Check template greater that another template.

        Parameters
        ----------
        t : Template
            Another template.

        Returns
        -------
        bool
            True - if template is greater than another template,
            False - otherwise.
        """

        return (self.lo, self.hi) > (t.lo, t.hi)

    #-----------------------------------------------------------------------------------------------

    def __le__(self, t):
        """
        Check template less or equal than another template.

        Parameters
        ----------
        t : Template
            Another template.

        Returns
        -------
        bool
            True - if template is less or equal than another template,
            False - otherwise.
        """

        return not (self > t)

    #-----------------------------------------------------------------------------------------------

    def __lt__(self, t):
        """
        Check template if less than another template.

        Parameters
        ----------
        t : Template
            Another template.

        Returns
        -------
        bool
            True - if template is less than another template,
            False - otherwise.
        """

        return (self.lo, self.hi) < (t.lo, t.hi)

    #-----------------------------------------------------------------------------------------------

    def __ge__(self, t):
        """
        Check template greater or equal than another template.

        Parameters
        ----------
        t : Template
            Another template.

        Returns
        -------
        bool
            True - is template is greater or equal than another template,
            False - otherwise.
        """

        return not (self < t)

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

    @staticmethod
    def create_random(n, step_lo, step_hi):

        # Create empty templates.
        ts = Templates()
        cur = 0

        # Adding templates
        for _ in range(n):
            off = random.randint(step_lo, step_hi)
            lo = cur + off

            # Create template.
            typ = random.randint(0, 3)
            if typ == 0:
                t = Template(lo, lo + 2, 1, -1)
            elif typ == 1:
                t = Template(lo, lo + 2, -1, -1)
            elif typ == 2:
                t = Template(lo, lo + 3, 3, -1)
            elif typ == 3:
                t = Template(lo, lo + 3, -3, -1)
            else:
                raise Exception('wrong type of generated template')

            ts.add(t)
            cur = t.lo

        return ts

    #-----------------------------------------------------------------------------------------------

    def __getitem__(self, i):
        """
        Get item.

        Parameters
        ----------
        i : int
            Index.

        Returns
        -------
        Template
            i-th template.
        """

        return self.items[i]

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
            if (self.items != []) and (t <= self.items[-1]):
                raise Exception('wrong order of adding templates')
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

    def is_empty(self):
        """
        Check if templates list is empty.

        Returns
        -------
        bool
            True - if templates list is empty,
            False - otherwise.
        """

        return self.count() == 0

    #-----------------------------------------------------------------------------------------------

    def border_length(self):
        """
        Length of border.

        Returns
        -------
        int
            Length of border.
        """

        if self.is_empty():
            return 0
        else:
            return self[-1].hi - self[0].lo

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
    ct = ts[ci]
    b[ci][0][0] = 0
    b[ci][ct.ubaln][1] = -1

    # Process all templates from last to first.
    while True:

        # Shift templates.
        ni, nt = ci, ct
        ci = ci - 1
        if ci < 0:
            break
        ct = ts[ci]

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

#---------------------------------------------------------------------------------------------------

def best_reduce_border_from_b_matrix(b):
    """
    Get best reduce border value from B matrix.

    Parameters
    ----------
    b : np.array
        3d-matrix.

    Returns
    -------
    int
        Best reduce border length value.
    """

    return min(b[0][0][0], b[0][0][1])

#---------------------------------------------------------------------------------------------------

def find_templates_from_b_matrix(b, ts):
    """
    Find templates from B matrix.

    Parameters
    ----------
    b : np.array
        3d matrix.
    ts : Templates
        Templates.

    Returns
    -------
    Templates
        Templates for use.
    """

    best = best_reduce_border_from_b_matrix(b)
    n = b.shape[0]

    # Horizontal position of route track.
    j = 0

    # Templates for return.
    ts_out = Templates()

    for i in range(n):
        if b[i][j][0] != best:
            t = ts[i]
            j = j - t.ubaln
            best = best + 1
            ts_out.add(t)

    return ts_out

#---------------------------------------------------------------------------------------------------

def smooth_border(ts):
    """
    Smooth border.

    Parameters
    ----------
    ts : Templates
        Templates.

    Returns
    -------
    float
        Border length reduce ratio.
    """

    b = build_b_matrix(ts)
    best = best_reduce_border_from_b_matrix(b)

    return -best / ts.border_length()

#===================================================================================================

def test():
    """
    Test.
    """

    # Test border smoothing.
    ts = Templates([Template(0, 2, 1, -1),
                    Template(1, 4, -3, -1),
                    Template(3, 5, 1, -1),
                    Template(5, 7, -1, -1),
                    Template(6, 8, 1, -1),
                    Template(7, 10, -3, -1),
                    Template(9, 11, 1, -1),
                    Template(10, 12, -1, -1),
                    Template(11, 14, +3, -1),
                    Template(13, 15, -1, -1),
                    Template(14, 16, +1, -1)])

    b = build_b_matrix(ts, False)
    best = best_reduce_border_from_b_matrix(b)

    assert best == -4

    ts_out = find_templates_from_b_matrix(b, ts)

    assert (ts_out[0].lo, ts_out[0].hi) == (5, 7)
    assert (ts_out[1].lo, ts_out[1].hi) == (7, 10)
    assert (ts_out[2].lo, ts_out[2].hi) == (11, 14)
    assert (ts_out[3].lo, ts_out[3].hi) == (14, 16)

#---------------------------------------------------------------------------------------------------

def test_random():
    """
    Random test.
    """

    max_step_lo = 10
    step_delta = 5
    templates_count = 100
    repeats_count = 10

    for step_lo in range(1, max_step_lo + 1):
        for step_hi in range(step_lo, step_lo + step_delta):
            r = [smooth_border(Templates.create_random(templates_count, step_lo, step_hi))
                 for _ in range(repeats_count)]
            print(f'steps = [{step_lo}, {step_hi}], r = {sum(r) / len(r)}')

#===================================================================================================

if __name__ == '__main__':
    test()

#===================================================================================================
