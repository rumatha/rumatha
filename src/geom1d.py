"""
1D geometry.
"""

#===================================================================================================

class Segment:
    """
    Segment.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, lo, hi):
        """
        Constructor.

        Parameters
        ----------
        lo : float
            Lo value.
        hi : float
            Hi value.
        """

        if lo > hi:
            raise Exception('geom1d:Segment.__init__: unable to create segment with lo > hi.')

        self.lo = lo
        self.hi = hi

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def whole():
        """
        Return segment of whole line.

        Returns
        -------
        Segment
            Segment [-inf, inf].
        """

        return Segment(float('-inf'), float('inf'))

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'[{self.lo} - {self.hi}]'

    #-----------------------------------------------------------------------------------------------

    @property
    def len(self):
        """
        Length.

        Returns
        -------
        float
            Length.
        """

        return self.hi - self.lo

    #-----------------------------------------------------------------------------------------------

    def mid(self):
        """
        Middle value.

        Returns
        -------
        float
            Middle value.
        """

        return 0.5 * (self.lo + self.hi)

    #-----------------------------------------------------------------------------------------------

    def split(self):
        """
        Split segment.

        Returns
        -------
        (Segment, Segment)
            Pair of segments.
        """

        return Segment(self.lo, self.mid()), Segment(self.mid(), self.hi)

#===================================================================================================

if __name__ == '__main__':
    pass

#===================================================================================================
