"""
3D geometry with real coordinated.
"""

import geom1d
import numpy as np

#===================================================================================================

class Point:
    """
    Point.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, x=0.0, y=0.0, z=0.0):
        """
        Constructor.

        Parameters
        ----------
        x : float
            X coordinate.
        y : float
            Y coordinate.
        z : float
            Z coordinate.
        """

        self.x = x
        self.y = y
        self.z = z

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def from_real_array(ar):
        """
        Create point from real array.

        Parameters
        ----------
        ar : np.array
            Array of coordinates.

        Returns
        -------
        Point
            Constructed point.
        """

        return Point(ar[0], ar[1], ar[2])

    #-----------------------------------------------------------------------------------------------

    def __sub__(self, p):
        """
        Subtraction of points.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        Vector
            Result vector.
        """

        return Vector(self.x - p.x, self.y - p.y, self.z - p.z)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def mid(p1, p2):
        """
        Middle point.

        Parameters
        ----------
        p1 : Point
            First point.
        p2 : Point
            Second point.

        Returns
        -------
        Point
            Middle point.
        """

        return Point(0.5 * (p1.x + p2.x), 0.5 * (p1.y + p2.y), 0.5 * (p1.z + p2.z))

#===================================================================================================

class Vector(Point):
    """
    Vector.
    """

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def dot(v1, v2):
        """
        Dot product.

        Parameters
        ----------
        v1 : Vector
            First vector.
        v2 : Vector
            Second vector.

        Returns
        -------
        float
            Result of dor product.
        """

        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z

#===================================================================================================

class Box:
    """
    Box.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, sx, sy, sz, eps):
        """
        Constructor.

        Parameters
        ----------
        sx : geom1d.Segment
            X segment.
        sy : geom1d.Segment
            Y segment.
        sz : geom1d.Segment
            Z segment.
        eps : float
            Epsilon.
        """

        self.sx = sx
        self.sy = sy
        self.sz = sz
        self.eps = eps

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def from_floats(xlo, xhi, ylo, yhi, zlo, zhi, eps):
        """
        Create from floats.

        Parameters
        ----------
        xlo : float
            X lo value.
        xhi : float
            Y hi value.
        ylo : float
            Y lo value.
        yhi : float
            Y hi value.
        zlo : float
            Z lo value.
        zhi : float
            Z hi value.
        eps : float
            Epsilon.

        Returns
        -------
        Box
            Created box.
        """

        return Box(geom1d.Segment(xlo, xhi),
                   geom1d.Segment(ylo, yhi),
                   geom1d.Segment(zlo, zhi), eps)

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'[{self.sx} x {self.sy} x {self.sz}, eps {self.eps}]'

    #-----------------------------------------------------------------------------------------------

    def mids(self):
        """
        Get mid values.

        Returns
        -------
        [float, float, float]
            Mids.
        """

        return [self.sx.mid(), self.sy.mid(), self.sz.mid()]

    #-----------------------------------------------------------------------------------------------

    def split(self, d):
        """
        Split box into two boxes.

        Parameters
        ----------
        d : str
            Direction ('x', 'y', 'z').

        Returns
        -------
        (Box, Box)
            Pair of boxes.
        """

        if d == 'x':
            sx1, sx2 = self.sx.split()
            return Box(sx1, self.sy, self.sz, self.eps), Box(sx2, self.sy, self.sz, self.eps)
        elif d == 'y':
            sy1, sy2 = self.sy.split()
            return Box(self.sx, sy1, self.sz, self.eps), Box(self.sx, sy2, self.sz, self.eps)
        else:
            if d != 'z':
                raise Exception(f'geom3d:Box.split: wrong direction of split ({d}).')
            sz1, sz2 = self.sz.split()
            return Box(self.sx, self.sy, sz1, self.eps), Box(self.sx, self.sy, sz2, self.eps)

#===================================================================================================

if __name__ == '__main__':
    pass

#===================================================================================================
