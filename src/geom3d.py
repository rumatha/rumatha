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

    def __init__(self, xlo, xhi, ylo, yhi, zlo, zhi):
        """
        Constructor.

        Parameters
        ----------
        xlo : float
            X lo value.
        xhi : float
            X hi value.
        ylo : float
            Y lo value.
        yhi : float
            Y hi value.
        zlo : float
            Z lo value.
        zhi : float
            Z hi value.
        """

        self.sx = geom1d.Segment(xlo, xhi)
        self.sy = geom1d.Segment(ylo, yhi)
        self.sz = geom1d.Segment(zlo, zhi)

#===================================================================================================

if __name__ == '__main__':
    pass

#===================================================================================================
