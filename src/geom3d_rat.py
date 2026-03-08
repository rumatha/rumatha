"""
3D geometry in rational coordinates.
"""

import numpy as np
from fractions import Fraction as Fr
import geom1d
import geom3d
import geom2d_rat
import bvh_tree
import matplotlib.pyplot as plt
import time

#===================================================================================================

class Point:
    """
    Point with three coordinates.
    """

    # Count of points created.
    counter = 0

    #-----------------------------------------------------------------------------------------------

    def __init__(self, x=Fr(0), y=Fr(0), z=Fr(0)):
        """
        Constructor.

        Parameters
        ----------
        x : Fraction
            x-coordinate.
        y : Fraction
            y-coordinate.
        z : Fraction
            z-coordinate.
        """

        # Increment counter.
        Point.counter = Point.counter + 1

        self.x = x
        self.y = y
        self.z = z

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def from_real_coords(x, y, z, denom):
        """
        Constructor from real coordinates.

        Parameters
        ----------
        x : float
            Real coordinate X.
        y : float
            Real coordinate Y.
        z : float
            Real coordinate Z.
        denom : int
            Denominator.

        Returns
        -------
        Point
            Constructed point.
        """

        return Point(Fr(int(x * denom), denom),
                     Fr(int(y * denom), denom),
                     Fr(int(z * denom), denom))

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def from_real_array(ar, denom):
        """
        Constructor from real array.

        Parameters
        ----------
        ar : np.array
            Array of coordinates.
        denom : int
            Denominator.

        Returns
        -------
        Point
            Constructed point.
        """

        return Point.from_real_coords(ar[0], ar[1], ar[2], denom)

    #-----------------------------------------------------------------------------------------------

    def get_real_array(self):
        """
        Get real array from rat coordinates.

        Returns
        -------
        np.array
            Array of coordinates.
        """

        coords = [self.x, self.y, self.z]

        return [f.numerator / f.denominator for f in coords]

    #-----------------------------------------------------------------------------------------------

    def __eq__(self, p):
        """
        Check equal with another point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if equal to another point,
            False - otherwise.
        """

        return (self.x == p.x) and (self.y == p.y) and (self.z == p.z)

    #-----------------------------------------------------------------------------------------------

    def __ne__(self, p):
        """
        Check not equal to another point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if not equal to another point.
            False - otherwise.
        """

        return not (self == p)

    #-----------------------------------------------------------------------------------------------

    def __ge__(self, p):
        """
        Check for GE with another point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if self >= p,
            False - otherwise.
        """

        if self.x > p.x:
            return True
        elif self.x < p.x:
            return False
        elif self.y > p.y:
            return True
        elif self.y < p.y:
            return False
        else:
            return self.z >= p.z

    #-----------------------------------------------------------------------------------------------

    def __gt__(self, p):
        """
        Check for GT with another point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if self > p,
            False - otherwise.
        """

        if self.x > p.x:
            return True
        elif self.x < p.x:
            return False
        elif self.y > p.y:
            return True
        elif self.y < p.y:
            return False
        else:
            return self.z > p.z

    #-----------------------------------------------------------------------------------------------

    def __le__(self, p):
        """
        Check for LE with another point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if self <= p,
            False - otherwise.
        """

        return not (self > p)

    #-----------------------------------------------------------------------------------------------

    def __lt__(self, p):
        """
        Check for LT with another point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if self < p,
            False - otherwise.
        """

        return not (self >= p)

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        x, y, z = self.x, self.y, self.z

        return f'P({x}={x.numerator/x.denominator}, '\
               f'{y}={y.numerator/y.denominator}, {z}={z.numerator/z.denominator})'

    #-----------------------------------------------------------------------------------------------

    def draw(self, plt, color='black', size=20):
        """
        Draw on plit.

        Parameters
        ----------
        plt : Plot
            Plot.
        color : str
            Color.
        size : int
            Size.
        """

        # Ignore no size points.
        if size > 0:

            # Ignore z coordinate.
            plt.scatter(self.x, self.y, color=color, s=size)

    #-----------------------------------------------------------------------------------------------

    def __sub__(self, p):
        """
        Two points difference.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        Vector
            Two points difference.
        """

        return Vector(self.x - p.x, self.y - p.y, self.z - p.z)

    #-----------------------------------------------------------------------------------------------

    def is_segment_end(self, s):
        """
        Check if point is end of segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if point is segment end,
            False - otherwise.
        """

        return (self == s.A) or (self == s.B)

    #-----------------------------------------------------------------------------------------------

    def is_triangle_vertex(self, t):
        """
        Check it point is triangle vertex.

        Parameters
        ----------
        t : Triangle
            Triangle.

        Returns
        -------
        bool
            True - if point is triangle vertex,
            False - otherwise.
        """

        return (self == t.A) or (self == t.B) or (self == t.C)

    #-----------------------------------------------------------------------------------------------

    def is_incident(self, obj):
        """
        Check if point is incident to obj (segment or triangle).

        Parameters
        ----------
        obj : Segment | Triangle
            Object.

        Returns
        -------
        bool
            True - if point is incident to object,
            False - otherwise.
        """

        if isinstance(obj, Segment):
            return self.is_segment_end(obj)
        else:
            assert isinstance(obj, Triangle)
            return self.is_triangle_vertex(obj)

    #-----------------------------------------------------------------------------------------------

    def projection_OXY(self):
        """
        Projection on OXY (ignore Z coord).

        Returns
        -------
        geom2d_rat.Point
            Projection.
        """

        return geom2d_rat.Point(self.x, self.y)

    #-----------------------------------------------------------------------------------------------

    def projection_OXZ(self):
        """
        Projection on OXZ (ignore Y coord).

        Returns
        -------
        geom2d_rat.Point
            Projection.
        """

        return geom2d_rat.Point(self.x, self.z)

    #-----------------------------------------------------------------------------------------------

    def projection_OYZ(self):
        """
        Projection of OYZ (ignore X coord).

        Returns
        -------
        geom2d_rat.Point
            Projection.
        """

        return geom2d_rat.Point(self.y, self.z)

#===================================================================================================

class Vector(Point):
    """
    Class vector.
    """

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'V({self.x}, {self.y}, {self.z})'

    #-----------------------------------------------------------------------------------------------

    def __add__(self, v):
        """
        Add two vectors.

        Parameters
        ----------
        v : Vector.
            Vector.

        Returns
        -------
        Vector
            New vector.
        """

        return Vector(self.x + v.x, self.y + v.y, self.z + v.z)

    #-----------------------------------------------------------------------------------------------

    def __sub__(self, v):
        """
        Subtract two vectors.

        Parameters
        ----------
        v : Vector
            Vector.

        Returns
        -------
        Vector
            New vector.
        """

        return Vector(self.x - v.x, self.y - v.y, self.z - v.z)

    #-----------------------------------------------------------------------------------------------

    def is_null(self):
        """
        Check if vector is null.

        Returns
        -------
        bool
            True - if vector is null,
            False - otherwise.
        """

        return (self.x == 0) and (self.y == 0) and (self.z == 0)

    #-----------------------------------------------------------------------------------------------

    def mod2(self):
        """
        Square of module.

        Returns
        -------
        Fraction
            Square of module.
        """

        return Vector.dot(self, self)

    #-----------------------------------------------------------------------------------------------

    def negate(self):
        """
        Negate vector.
        """

        self.x = -self.x
        self.y = -self.y
        self.z = -self.z

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
        Fraction
            Result.
        """

        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def vector_product(a, b):
        """
        Vector product.

        Parameters
        ----------
        a : Vector
            First vector.
        b : Vector
            Second vector.

        Returns
        -------
        Vector
            Result vector.
        """

        return Vector(a.y * b.z - a.z * b.y,
                      a.z * b.x - a.x * b.z,
                      a.x * b.y - a.y * b.x)

#===================================================================================================

class Points:
    """
    Points.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        self.items = []

    #-----------------------------------------------------------------------------------------------

    def __getitem__(self, i):
        """
        Get i-th point.

        Parameters
        ----------
        i : int
            Index.

        Returns
        -------
        Point
            Point.
        """

        return self.items[i]

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Points[{self.items}]'

    #-----------------------------------------------------------------------------------------------

    def draw(self, plt, color='black', size=20):
        """
        Draw on plot.

        Parameters
        ----------
        plt : Plot
            Plot.
        color : str
            Color.
        size : int
            Size.
        """

        for p in self.items:
            p.draw(plt, color=color, size=size)

    #-----------------------------------------------------------------------------------------------

    def count(self):
        """
        Count of points.

        Returns
        -------
        int
            Count of points.
        """

        return len(self.items)

    #-----------------------------------------------------------------------------------------------

    def add(self, p):
        """
        Add new point.

        Parameters
        ----------
        p : Point
            Point.
        """

        self.items.append(p)

    #-----------------------------------------------------------------------------------------------

    def is_contain_point(self, p):
        """
        Check if points set contains point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if point in points set,
            False - otherwise.
        """

        for pi in self.items:
            if pi == p:
                return True

        return False

    #-----------------------------------------------------------------------------------------------

    def add_unique(self, p):
        """
        Add new unique point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if point was added,
            False - if point was not added.
        """

        if self.is_contain_point(p):
            return False

        self.add(p)

        return True

    #-----------------------------------------------------------------------------------------------

    def sort(self):
        """
        Sort points.
        """

        self.items.sort()

#===================================================================================================

class Line:
    """
    Line in space.
    """

    # Counter for created lines.
    counter = 0

    #-----------------------------------------------------------------------------------------------

    def __init__(self, x0, y0, z0, m, n, p):
        """
        Line constructor.

        x = x0 + tm
        y = y0 + tn
        z = z0 + tp

        Parameters
        ----------
        x0 : Fraction
            X coordinate of base point.
        y0 : Fraction
            Y coordinate of base point.
        z0 : Fraction
            Z coordinate of base point.
        m : Fraction
            Parameter for X direction.
        n : Fraction
            Parameter for Y direction.
        p : Fraction
            Parameter for Z direction.
        """

        # Normalize vector.
        if m != 0:
            n = n / m
            p = p / m
            m = Fr(1)
        elif n != 0:
            p = p / n
            n = Fr(1)
        else:
            if p == 0:
                raise Exception('geom3d_rat:Line.__init__: constructor from zero vector m, n, p.')
            p = Fr(1)

        # Normalize point.
        # Normalization is needed for fast comparing of two lines.
        if m != 0:
            t = -x0 / m
            y0 = y0 + t * n
            z0 = z0 + t * p
            x0 = Fr(0)
        elif n != 0:
            t = -y0 / n
            x0 = x0 + t * m
            z0 = z0 + t * p
            y0 = Fr(0)
        else:
            t = -z0 / p
            x0 = x0 + t * m
            y0 = y0 + t * n
            z0 = Fr(0)

        self.P0 = Point(x0, y0, z0)
        self.v = Vector(m, n, p)

        # Increment counter.
        Line.counter = Line.counter + 1

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def from_point_and_vector(p, v):
        """
        Constructor from point and vector.

        Parameters
        ----------
        p : Point
            Base point.
        v : Vector
            Direction vector.

        Returns
        -------
        Line
            Constructed line.
        """

        return Line(p.x, p.y, p.z, v.x, v.y, v.z)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def from_points(p1, p2):
        """
        Constructor from two points.

        Parameters
        ----------
        p1 : Point
            First point.
        p2 : Point
            Second point.

        Returns
        -------
        Line
            Constructed line.
        """

        return Line.from_point_and_vector(p1, p2 - p1)

    #-----------------------------------------------------------------------------------------------

    def __eq__(self, line):
        """
        Check equal.

        Parameters
        ----------
        line : Line
            Line.

        Returns
        -------
        bool
            True - if lines are equal,
            False - otherwise.
        """

        return (self.P0 == line.P0) and (self.v == line.v)

    #-----------------------------------------------------------------------------------------------

    def __ne__(self, line):
        """
        Check not equal.

        Parameters
        ----------
        line : Line
            Line.

        Returns
        -------
        bool
            True - if lines are not equal,
            False - otherwise.
        """

        return not (self == line)

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Line({self.P0} + t * {self.v})'

    #-----------------------------------------------------------------------------------------------

    def is_have_point(self, p):
        """
        Check if line has point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if line has point,
            False - otherwise.
        """

        dx, dy, dz = p.x - self.P0.x, p.y - self.P0.y, p.z - self.P0.z
        m, n, p = self.v.x, self.v.y, self.v.z

        # Check (dx, dy, dz) codirected with (m, n, p)

        # Find coefficient.
        if m != 0:
            t = dx / m
        elif n != 0:
            t = dy / n
        else:
            if p == 0:
                raise Exception('geom3d_rat:Line.is_have_point: incorrent line vector.')
            t = dz / p

        return (m * t == dx) and (n * t == dy) and (p * t == dz)

#===================================================================================================

class Segment:
    """
    Segment in space.
    """

    # Counter for created segments.
    counter = 0

    #-----------------------------------------------------------------------------------------------

    def __init__(self, A, B):
        """
        Constructor by points.

        Parameters
        ----------
        A : Point
            A Point.
        B : Point
            B Point.
        """

        # Construction of zero length segment is forbidden.
        if A == B:
            raise Exception('geom3d_rat:Segment:__init__: construct segments of zero length.')

        self.A = A
        self.B = B

        # Keep points in sorted way.
        self.sort_points()

        # Construct line for this segment.
        self.line = Line.from_points(self.A, self.B)

        # Increment counter.
        Segment.counter = Segment.counter + 1

    #-----------------------------------------------------------------------------------------------

    def __eq__(self, s):
        """
        Check equal to another segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if equal to another segment,
            False - otherwise.
        """

        # Points are sorted.
        return (self.A == s.A) and (self.B == s.B)

    #-----------------------------------------------------------------------------------------------

    def __ne__(self, s):
        """
        Check not equal to another segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if not equal to another segment,
            False - otherwise.
        """

        return not (self == s)

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Segm[{self.A}, {self.B}]'

    #-----------------------------------------------------------------------------------------------

    def draw(self, plt, color='black', linewidth='2', size=20):
        """
        Draw on plot.

        Parameters
        ----------
        plt : Plot
            Plot.
        color : str
            Color.
        linewidth : str
            Line width.
        size : int
            Size.
        """

        # We ignore z coordinate, draw only on OXY plane.
        x = [self.A.x, self.B.x]
        y = [self.A.y, self.B.y]

        # Do not draw segment with zero width.
        if linewidth != '0':
            plt.plot(x, y, color=color, linewidth=linewidth)

        # Draw ends.
        self.A.draw(plt, color=color, size=size)
        self.B.draw(plt, color=color, size=size)

    #-----------------------------------------------------------------------------------------------

    def mod2(self):
        """
        Square of module.

        Returns
        -------
        Fraction
            Square of module.
        """

        return (self.A.x - self.B.x)**2 + (self.A.y - self.B.y)**2 + (self.A.z - self.B.z)**2

    #-----------------------------------------------------------------------------------------------

    def is_triangle_side(self, t):
        """
        Check if segment is triangle side.

        Parameters
        ----------
        t : Triangle
            Triangle.

        Returns
        -------
        bool
            True - if segment is triangle side,
            False - otherwise.
        """

        return (self == t.AB) or (self == t.BC) or (self == t.AC)

    #-----------------------------------------------------------------------------------------------

    def is_incident(self, obj):
        """
        Check if segment is incident to object (point or triangle).

        Parameters
        ----------
        obj : Point | Triangle
            Object.

        Returns
        -------
        bool
            True - if segment is incident to point or triangle,
            False - otherwise.
        """

        if isinstance(obj, Point):
            return obj.is_segment_end(self)
        else:
            assert isinstance(obj, Triangle)
            return self.is_triangle_side(obj)

    #-----------------------------------------------------------------------------------------------

    def is_point_between_ends(self, p):
        """
        Check if point is between ends.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if point is between ends.
            False - otherwise.
        """

        return (self.A <= p) and (p <= self.B)

    #-----------------------------------------------------------------------------------------------

    def is_have_point(self, p):
        """
        Check if segment has point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if segment has point.
            False - otherwise.
        """

        return self.line.is_have_point(p) and self.is_point_between_ends(p)

    #-----------------------------------------------------------------------------------------------

    def is_adjacent(self, s):
        """
        Check for adjacency with another segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if adjacent with another segment,
            False - otherwise.
        """

        # Same segments are not adjacent.
        if self == s:
            return False

        # Check for equal lines.
        if self.line == s.line:
            # Segments are adjacent if the second end of one segment is the first end of another.
            return (self.B == s.A) or (s.B == self.A)

        # Two segments of general placement are adjacent if they have common end.
        return self.A.is_segment_end(s) or self.B.is_segment_end(s)

    #-----------------------------------------------------------------------------------------------

    def is_conflict(self, s):
        """
        Check if segment conflicts with another segment.
        Segments are called conflicted if they have common point
        which is inner point of one or both segments.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if segment conflicts with another segment,
            False - otherwise.
        """

        # Adjacent segments do not conflict.
        if self.is_adjacent(s):
            return False

        # Find intersection.
        r = Intersection.segment_segment(self, s)

        # No intersection - no conflict.
        return not (r is None)

    #-----------------------------------------------------------------------------------------------

    def sort_points(self):
        """
        Sort points.
        """

        if self.A > self.B:
            self.A, self.B = self.B, self.A

    #-----------------------------------------------------------------------------------------------

    def split(self, ps):
        """
        Split segment into list of segments.

        Parameters
        ----------
        ps : Points
            Points.

        Returns
        -------
        [Segment]
            List of segments.
        """

        ps.sort()

        i, n = 0, ps.count()

        # Find first point inside segment.
        while (i < n) and (ps[i] <= self.A):
            i = i + 1

        # If all points are not greater than A then terminate.
        if i == n:
            return self

        # if first point greater than A also greater than B, then terminate.
        if ps[i] >= self.B:
            return self

        # There is points inside segment.
        j = n - 1
        while ps[j] >= self.B:
            j = j - 1

        # Separate points are ps[i] < .. < ps[j].
        ss = []
        prev = self.A
        while i <= j:
            ss.append(Segment(prev, ps[i]))
            if i == j:
                ss.append(Segment(ps[i], self.B))
                break
            prev = ps[i]
            i = i + 1
        return ss

#===================================================================================================

class Segments:
    """
    Segments.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        self.items = []

    #-----------------------------------------------------------------------------------------------

    def __getitem__(self, i):
        """
        Get i-th segment.

        Parameters
        ----------
        i : int
            Index.

        Returns
        -------
        Segment
            Segment on i-th position.
        """

        return self.items[i]

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Segments[{self.items}]'

    #-----------------------------------------------------------------------------------------------

    def draw(self, plt, color='black', linewidth='2', size=20):
        """
        Draw on plot.

        Parameters
        ----------
        plt : Plot
            Plot.
        color : str
            Color.
        linewidth : str
            Line width.
        size : int
            Size.
        """

        for s in self.items:
            s.draw(plt, color=color, linewidth=linewidth, size=size)

    #-----------------------------------------------------------------------------------------------

    def count(self):
        """
        Count of segments.

        Returns
        -------
        int
            Count of segments.
        """

        return len(self.items)

    #-----------------------------------------------------------------------------------------------

    def add(self, s):
        """
        Add segment.

        Parameters
        ----------
        s : Segment
            Segment.
        """

        self.items.append(s)

    #-----------------------------------------------------------------------------------------------

    def is_contain_segment(self, s):
        """
        Check if segments set contains segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if segment is in segments set,
            False - otherwise.
        """

        for si in self.items:
            if si == s:
                return True

        return False

    #-----------------------------------------------------------------------------------------------

    def is_contain_point_as_segment_end(self, p):
        """
        Check if point is some segment end.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if point is end of some segment,
            False - otherwise.
        """

        for si in self.items:
            if p.is_segment_end(si):
                return True

        return False

    #-----------------------------------------------------------------------------------------------

    def add_unique(self, s):
        """
        Add unique segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if new segment was added,
            False - if segment was not added.
        """

        if self.is_contain_segment(s):
            return False

        self.add(s)

        return True

    #-----------------------------------------------------------------------------------------------

    def adds(self, ss):
        """
        Add list of segments.

        Parameters
        ----------
        ss : [Segment]
            List of segments.
        """

        for s in ss:
            self.add(s)

    #-----------------------------------------------------------------------------------------------

    def adds_unique(self, ss):
        """
        Add list of segments as unique segments.

        Parameters
        ----------
        ss : [Segment]
            List of segments.
        """

        for s in ss:
            self.add_unique(s)

    #-----------------------------------------------------------------------------------------------

    def points(self):
        """
        Get all points.

        Returns
        -------
        Points
            Points.
        """

        ps = Points()

        for s in self.items:
            ps.add_unique(s.A)
            ps.add_unique(s.B)

        return ps

    #-----------------------------------------------------------------------------------------------

    def sort(self, fun):
        """
        Sort set of segments.

        Parameters
        ----------
        fun : fun
            Key function.
        """

        self.items.sort(key=fun)

    #-----------------------------------------------------------------------------------------------

    def is_conflict(self, s):
        """
        Check if segments conflict with given segment.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if segments conflict with given segment,
            False - otherwise.
        """

        for si in self.items:
            if si.is_conflict(s):
                return True

        return False

    #-----------------------------------------------------------------------------------------------

    def is_have_conflict(self):
        """
        Check if there is segment-segment conflict.

        Returns
        -------
        bool
            True - if there is segment-segment conflict,
            False - otherwise.
        """

        n = self.count()

        for i in range(n):
            for j in range(i + 1, n):
                if self[i].is_conflict(self[j]):
                    return True

        return False

    #-----------------------------------------------------------------------------------------------

    def split_segment_by_point(self, i, p):
        """
        Split segment by point.

        Parameters
        ----------
        i : int
            Segment position.
        p : Point
            Point.
        """

        s = self[i]

        # If point is segment end then nothing to do.
        if p.is_segment_end(s):
            return

        # Split segment into two parts.
        # First part substitutes initial segment, second part goes to the end of segments list.
        s1, s2 = Segment(s.A, p), Segment(p, s.B)
        self[i] = s1
        # Just add it because segment s2 can cause new conflicts.
        self.add(s2)

    #-----------------------------------------------------------------------------------------------

    def fix_conflicts(self):
        """
        Fix conflicts.
        """

        # Create segments and points for cutting.
        ss = [(s, Points()) for s in self.items]
        self.items = []

        # Check all pairs of segments.
        n = len(ss)
        for i in range(n):
            (s1, int1) = ss[i]
            for j in range(i + 1, n):
                (s2, int2) = ss[j]
                r = Intersection.segment_segment(s1, s2)
                if isinstance(r, Point):
                    int1.add_unique(r)
                    int2.add_unique(r)
                elif isinstance(r, Segment):
                    int1.add_unique(r.A)
                    int1.add_unique(r.B)
                    int2.add_unique(r.A)
                    int2.add_unique(r.B)

        # Split all and add back.
        for (s, int) in ss:
            parts = s.split(int)
            if isinstance(parts, Segment):
                self.add(parts)
            else:
                self.adds(parts)

    #-----------------------------------------------------------------------------------------------

    def triangles(self):
        """
        Get triangles.

        Returns
        -------
        Triangles
            Get triangles.
        """

        ts = Triangles()
        n = self.count()

        # Check all triangles.
        for i in range(n):
            s = self[i]
            for j in range(i + 1, n):
                q = self[j]
                if s.is_adjacent(q):
                    for k in range(j + 1, n):
                        r = self[k]
                        if r.is_adjacent(s):
                            if r.is_adjacent(q):
                                ps = Points()
                                for p in [s.A, s.B, q.A, q.B, r.A, r.B]:
                                    ps.add_unique(p)
                                if ps.count() == 3:
                                    t = Triangle(ps[0], ps[1], ps[2])
                                    ts.add(t)

        return ts

#===================================================================================================

class PointsAndSegments:
    """
    Class that holds points and segments.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        self.points = Points()
        self.segments = Segments()

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'P&S: {self.points} + {self.segments}'

    #-----------------------------------------------------------------------------------------------

    def draw(self, plt, color='black', linewidth='2', size=20):
        """
        Draw on plot.

        Parameters
        ----------
        plt : Plot
            Plot.
        color : str
            Color.
        linewidth : str
            Line width.
        size : int
            Size.
        """

        self.segments.draw(plt, color, linewidth, size)
        self.points.draw(plt, color, size)

    #-----------------------------------------------------------------------------------------------

    def points_count(self):
        """
        Points count.

        Returns
        -------
        int
            Points count.
        """

        return self.points.count()

    #-----------------------------------------------------------------------------------------------

    def segments_count(self):
        """
        Segments count.

        Returns
        -------
        int
            Segments count.
        """

        return self.segments.count()

    #-----------------------------------------------------------------------------------------------

    def is_empty(self):
        """
        Check if points and segments empty.

        Returns
        -------
        bool
            True - if set is empty,
            False - otherwise.
        """

        return (self.points_count() == 0) and (self.segments_count() == 0)

    #-----------------------------------------------------------------------------------------------

    def add_point(self, p):
        """
        Add point.

        Parameters
        ----------
        p : Point
            Point.
        """

        self.points.add(p)

    #-----------------------------------------------------------------------------------------------

    def add_segment(self, s):
        """
        Add segment.

        Parameters
        ----------
        s : Segment
            Segment.
        """

        self.segments.add(s)

    #-----------------------------------------------------------------------------------------------

    def add_unique_point(self, p):
        """
        Add unique point.

        We add point only if it is not in points set,
        and if it is not end of some segment.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if point was added (because it was unique),
            False - otherwise.
        """

        # Try to find point as object or segment end.
        if self.points.is_contain_point(p) or self.segments.is_contain_point_as_segment_end(p):
            return False

        self.add_point(p)

        return True

    #-----------------------------------------------------------------------------------------------

    def add_unique_segment(self, s):
        """
        Add new unique segment.

        We can add segment only if there is no such segment.
        If any point is end of new segment then delete this point.

        Parameters
        ----------
        s : Segment
            Segment.

        Returns
        -------
        bool
            True - if segment is added,
            False - otherwise.
        """

        # If there is such segment then do not add.
        if self.segments.is_contain_segment(s):
            return False

        # Check all points and write them into new set.
        ps = Points()
        for pi in self.points.items:
            if not pi.is_segment_end(s):
                ps.add(pi)
        self.points = ps

        self.segments.add(s)

        return True

    #-----------------------------------------------------------------------------------------------

    def adds_segments(self, ss):
        """
        Add segments.

        Parameters
        ----------
        ss : [Segment]
            Segments.
        """

        for s in ss:
            self.add(s)

    #-----------------------------------------------------------------------------------------------

    def adds_unique_segments(self, ss):
        """
        Add unique segments.

        Parameters
        ----------
        ss : [Segment]
            Segments.
        """

        for s in ss:
            self.add_unique_segment(s)

    #-----------------------------------------------------------------------------------------------

    def is_have_segment_point_conflict(self):
        """
        Check if there is segment-point conflict.
        Segment-point conflict is equal segment has point.

        Returns
        -------
        bool
            True - if there is segment-point conflict,
            False - otherwise.
        """

        for si in self.segments.items:
            for pi in self.points.items:
                if si.is_have_point(pi):
                    return True

        return False

    #-----------------------------------------------------------------------------------------------

    def is_have_segment_segment_conflict(self):
        """
        Check if there is segment-segment conflict.

        Returns
        -------
        bool
            True - if there is segment-segment conflict,
            False - otherwise.
        """

        return self.segments.is_have_conflict()

    #-----------------------------------------------------------------------------------------------

    def is_have_conflict(self):
        """
        Check if there is conflict.

        Returns
        -------
        bool
            True - if there is conflict,
            False - otherwise.
        """

        return self.is_have_segment_point_conflict() or self.is_have_segment_segment_conflict()

    #-----------------------------------------------------------------------------------------------

    def fix_conflicts(self):
        """
        Fix all conflicts.
        """

        self.segments.fix_conflicts()

    #-----------------------------------------------------------------------------------------------

    def points(self):
        """
        Get all points.
        Points and segments ends.

        Returns
        -------
        Points
            Points and segment ends.
        """

        ps = self.segments.points()

        for p in self.points.items:
            ps.add_unique(p)

        return ps

    #-----------------------------------------------------------------------------------------------

    def possible_segments(self):
        """
        Get all possible segments.
        Possible segment is constructed from pair of points.

        Returns
        -------
        Segments
            Possibble segments.
        """

        # Get all points.
        ps = self.segments.points()
        for p in self.points:
            ps.add_unique(p)
        n = ps.count()
        ss = Segments()

        # Construct segments.

        for i in range(n):
            for j in range(i + 1, n):
                ss.add(Segment(ps[i], ps[j]))

        return ss

    #-----------------------------------------------------------------------------------------------

    def minimal_segments_coverage(self):
        """
        Get minimal segments coverage.

        Returns
        -------
        Segments
            Minimal segments coverage.
        """

        # Get all possible segments.
        ps = self.possible_segments()

        # Create new segments container.
        ms = Segments()
        for s in self.segments:
            ms.add(s)

        # First get segments with small length.
        ps.sort(fun=lambda s: s.mod2())

        # Try to take as many segments as possible.
        # We do not need add unique segment because there is no conflict.
        for s in ps.items:
            if not ms.is_conflict(s):
                ms.add(s)

        return ms

#===================================================================================================

class Plane:
    """
    Plane.
    """

    # Counter of created planes.
    counter = 0

    #-----------------------------------------------------------------------------------------------

    def __init__(self, a, b, c, d):
        """
        Constructor.

        Parameters
        ----------
        a : Fraction
            a coefficient.
        b : Fraction
            b coefficient.
        c : Fraction
            c coefficient.
        d : Fraction
            d coefficient.
        """

        self.a = a
        self.b = b
        self.c = c
        self.d = d

        # Construct normal.
        self.normal = Vector(self.a, self.b, self.c)

        # Increment counter.
        Plane.counter = Plane.counter + 1

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def from_points(A, B, C):
        """
        Create plane from three points.

        Plane types.
        1) if a != 0 then a = 1
           x + by + cz + d = 0
        2) if a = 0 and b != 0 then b = 1
           b + cz + d = 0
        3) if a = 0 and b = 0 and c != 0 then c = 1
           z + d = 0
        4) if a = 0 and b = 0 and c = 0 then there is no plane.

        Parameters
        ----------
        A : Point
            A Point.
        B : Point
            B Point.
        C : Point
            C Point.

        Returns
        -------
        Plane
            Result plane.
        """

        x1, y1, z1, x2, y2, z2, x3, y3, z3 = A.x, A.y, A.z, B.x, B.y, B.z, C.x, C.y, C.z

        # src: https://guimc.bmstu.ru/wp-content/uploads/2018/11/lecture_2.1.pdf
        # vectors M1M  = (x - x1, y - y1, z - z1)
        #         M1M2 = (x2 - x1, y2 - y1, z2 - z1)
        #         M1M3 = (x3 - x1, y3 - y1, z3 - z1)
        # These vectors lie in one plane when they are complanar:
        #   | x - x1    y - y1    z - z1  |
        #   | x2 - x1   y2 - y1   z2 - z1 | = 0
        #   | x3 - x1   y3 - y1   z3 - z1 |
        x21, y21, z21 = x2 - x1, y2 - y1, z2 - z1
        x31, y31, z31 = x3 - x1, y3 - y1, z3 - z1

        #   | x - x1    y - y1    z - z1 |
        #   |  x21       y21       z21   | = 0
        #   |  x31       y31       z31   |
        # Calculate determinant:
        # (x - x1) * |y21 z21| - (y - y1) * |x21 z21| + (z - z1) * |x21 y21| = 0
        #            |y31 z31|              |x31 z31|              |x31 y31|
        dx = y21 * z31 - y31 * z21
        dy = -(x21 * z31 - x31 * z21)
        dz = (x21 * y31 - x31 * y21)

        # (x - x1) * dx + (y - y1) * dy + (z - z1) * dz = 0
        # x dx - x1 dx + y dy - y1 dy + z dz - z1 dz = 0
        # x dx + y dy + z dz + (-(x1 dx + y1 dy + z1 dz)) = 0
        a = dx
        b = dy
        c = dz
        d = -(x1 * dx + y1 * dy + z1 * dz)

        if a != 0:
            b = b / a
            c = c / a
            d = d / a
            a = Fr(1)
        elif b != 0:
            c = c / b
            d = d / b
            b = Fr(1)
        elif c != 0:
            d = d / c
            c = Fr(1)
        else:
            raise Exception(f'Plane can not be constructed from points {A}, {B}, {C}.')

        return Plane(a, b, c, d)

    #-----------------------------------------------------------------------------------------------

    def __eq__(self, p):
        """
        Check equal with another plane.

        Parameters
        ----------
        p : Plane
            Plane.

        Returns
        -------
        bool
            True - if equal to another plane,
            False - otherwise.
        """

        return (self.a == p.a) and (self.b == p.b) and (self.c == p.c) and (self.d == p.d)

    #-----------------------------------------------------------------------------------------------

    def __ne__(self, p):
        """
        Check not equal with another plane.

        Parameters
        ----------
        p : Plane
            Plane.

        Returns
        -------
        bool
            True - if not equal to another plane,
            False - otherwise.
        """

        return not (self == p)

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Plane({self.a} x + {self.b} y + {self.c} z + {self.d})'

    #-----------------------------------------------------------------------------------------------

    def val(self, p):
        """
        Value of point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        Fraction
            Value.
        """

        return self.a * p.x + self.b * p.y + self.c * p.z + self.d

    #-----------------------------------------------------------------------------------------------

    def is_have_point(self, p):
        """
        Check if plane has point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if plane has point,
            False - otherwise.
        """

        return self.val(p) == 0

    #-----------------------------------------------------------------------------------------------

    def is_two_points_strong_on_one_side(self, p1, p2):
        """
        Check if two points strong on one side.

        Parameters
        ----------
        p1 : Point
            First point.
        p2 : Point
            Second point.

        Returns
        -------
        bool
            True - if two points on one side of plane,
            False - otherwise.
        """

        return self.val(p1) * self.val(p2) > 0

    #-----------------------------------------------------------------------------------------------

    def is_three_points_strong_on_one_side(self, p1, p2, p3):
        """
        Check if three points strong on one side.

        Parameters
        ----------
        p1 : Point
            First point.
        p2 : Point
            Second point.
        p3 : Point
            Three point.

        Returns
        -------
        bool
            True - if three points on one side of plane,
            False - otherwise.

        """

        v1, v2, v3 = self.val(p1), self.val(p2), self.val(p3)

        return (v1 * v2 > 0) and (v1 * v3 > 0)

    #-----------------------------------------------------------------------------------------------

    def is_perpendicular_with_plane(self, pl):
        """
        Check if plane is perpendicular with plane.

        Parameters
        ----------
        pl : Plane
            Plane.

        Returns
        -------
        bool
            True - if perpendicular with plane,
            False - otherwise.
        """

        return Vector.dot(self.normal, pl.normal) == 0

#===================================================================================================

class Triangle:
    """
    Triangle in space.
    """

    # Counter for created triangles.
    counter = 0

    #-----------------------------------------------------------------------------------------------

    def __init__(self, A, B, C):
        """
        Constructor.

        Parameters
        ----------
        A : Point
            A Point.
        B : Point
            B Point.
        C : Point
            C Point.
        """

        if (A == B) or (B == C) or (A == C):
            raise Exception('geom3d_rat:Triangle.__init__: create triangle from points, '
                            f'which have duplicates ({A}, {B}, {C})')

        self.A = A
        self.B = B
        self.C = C
        self.points = [self.A, self.B, self.C]

        # Init sides at a time.
        self.AB = Segment(self.A, self.B)
        self.BC = Segment(self.B, self.C)
        self.AC = Segment(self.A, self.C)
        self.sides = [self.AB, self.BC, self.AC]

        # Plane.
        self.plane = Plane.from_points(self.A, self.B, self.C)

        # Construct outer normal for triangle.
        self.outer_normal = Vector.vector_product(self.B - self.A, self.C - self.A)

        # Increment counter.
        Triangle.counter = Triangle.counter + 1

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Tri<{self.A}, {self.B}, {self.C}>'

    #-----------------------------------------------------------------------------------------------

    def draw(self, plt, color='black', linewidth='2', size=20):
        """
        Draw on plot.

        Parameters
        ----------
        plt : Plot
            Plot.
        color : str
            Color.
        linewidth : str
            Color.
        size : int
            Size.
        """

        # First draw sides without points.
        for s in self.sides:
            s.draw(plt, color=color, linewidth=linewidth, size=0)

        # Then draw points over sides.
        for p in self.points:
            p.draw(plt, color=color, size=size)

    #-----------------------------------------------------------------------------------------------

    def is_incident(self, obj):
        """
        Check if triangle is incident to point or segment.

        Parameters
        ----------
        obj : Point | Segment
            Object.

        Returns
        -------
        bool
            True - if triangle is incident to point or segment,
            False - otherwise.
        """

        if isinstance(obj, Point):
            return obj.is_triangle_vertex(self)
        else:
            assert isinstance(obj, Segment)
            return obj.is_triangle_side(self)

    #-----------------------------------------------------------------------------------------------

    def is_perpendicular_with_plane(self, pl):
        """
        Check if perpendicular with plane.

        Parameters
        ----------
        pl : Plane
            Plane.

        Returns
        -------
        bool
            True - if perpendicular with plane,
            False - otherwise.
        """

        return self.plane.is_perpendicular_with_plane(pl)

    #-----------------------------------------------------------------------------------------------

    def is_conflict(self, t):
        """
        Check conflict with another triangle.

        Parameters
        ----------
        t : Triangle
            Triangle.

        Returns
        -------
        bool
            True - if there is conflict with another triangle,
            False - otherwise.
        """

        r = Intersection.triangle_triangle(self, t)

        return isinstance(r, Segments)

    #-----------------------------------------------------------------------------------------------

    def area2(self):
        """
        Squared area.

        Returns
        -------
        Fraction
            Squared area.
        """

        return Vector.vector_product(self.B - self.A, self.C - self.A).mod2() * Fr(1, 4)

    #-----------------------------------------------------------------------------------------------

    def flip_normal(self):
        """
        Flip normal.
        """

        self.B, self.C = self.C, self.B
        self.outer_normal.negate()

    #-----------------------------------------------------------------------------------------------

    def projection_OXY(self):
        """
        Projection on OXY (ignore Z coordinate).

        Returns
        -------
        geom2d_rat.Triangle
            Projection.
        """

        return geom2d_rat.Triangle(self.A.projection_OXY(),
                                   self.B.projection_OXY(),
                                   self.C.projection_OXY())

    #-----------------------------------------------------------------------------------------------

    def projection_OXZ(self):
        """
        Projection on OXZ (ignore Y coordinate).

        Returns
        -------
        geom2d_rat.Triangle
            Projection.
        """

        return geom2d_rat.Triangle(self.A.projection_OXZ(),
                                   self.B.projection_OXZ(),
                                   self.C.projection_OXZ())

    #-----------------------------------------------------------------------------------------------

    def projection_OYZ(self):
        """
        Projection on OYZ (ignore X coordinate).

        Returns
        -------
        geom2d_rat.Triangle
            Projection.
        """

        return geom2d_rat.Triangle(self.A.projection_OYZ(),
                                   self.B.projection_OYZ(),
                                   self.C.projection_OYZ())

    #-----------------------------------------------------------------------------------------------

    def is_have_point(self, p):
        """
        Check if triangle has point.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if triangle has point,
            False - otherwise.
        """

        if not self.is_perpendicular_with_plane(OXY):
            return self.projection_OXY().is_have_point(p.projection_OXY())
        elif not self.is_perpendicular_with_plane(OXZ):
            return self.projection_OXZ().is_have_point(p.projection_OXZ())
        elif not self.is_perpendicular_with_plane(OYZ):
            return self.projection_OYZ().is_have_point(p.projection_OYZ())
        else:
            assert False

    #-----------------------------------------------------------------------------------------------

    def is_have_point_on_sides(self, p):
        """
        Check if triangle has point on one of its sides.

        Parameters
        ----------
        p : Point
            Point.

        Returns
        -------
        bool
            True - if point is on side,
            False - otherwise.
        """

        return self.AB.is_have_point(p) or self.BC.is_have_point(p) or self.AC.is_have_point(p)

    #-----------------------------------------------------------------------------------------------

    def triangulate(self, pas):
        """
        Triangulate by set of points and segments.

        Parameters
        ----------
        pas : PointsAndSegments
            Set of points and segments.

        Returns
        -------
        Triangles
            Triangles.
        """

        # For triangulation we have to add triangle sides into points and segments set.
        # After this fix conflicts.
        pas.adds_unique_segments(self.sides)
        pas.fix_conflicts()

        if pas.is_have_conflict():
            raise Exception('Triangle.triangulate: points and segments conflict')

        # Get minimal segments coverage.
        msc = pas.minimal_segments_coverage()

        # Get all possible triangles.
        ts = msc.triangles()

        # Make mosaic from triangles.
        ts.construct_mosaic()

        return ts

#===================================================================================================

class Triangles:
    """
    Triangles class.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        self.items = []

    #-----------------------------------------------------------------------------------------------

    def __getitem__(self, i):
        """
        Get i-th triangle.

        Parameters
        ----------
        i : int
            Index.

        Returns
        -------
        Triangle
            Get i-th triangle.
        """

        return self.items[i]

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Triangles[{self.items}]'

    #-----------------------------------------------------------------------------------------------

    def draw(self, plt, color='black', linewidth='2', size=20):
        """
        Draw on plot.

        Parameters
        ----------
        plt : Plot
            Plot.
        color : str
            Color.
        linewidth : str
            Line width.
        size : int
            Point size.
        """

        for t in self.items:
            t.draw(plt, color=color, linewidth=linewidth, size=size)

    #-----------------------------------------------------------------------------------------------

    def count(self):
        """
        Count of triangles.

        Returns
        -------
        int
            Count of triangles.
        """

        return len(self.items)

    #-----------------------------------------------------------------------------------------------

    def add(self, t):
        """
        Add triangle.

        Parameters
        ----------
        t : Triangle
            Triangle.
        """

        self.items.append(t)

    #-----------------------------------------------------------------------------------------------

    def sort(self, fun):
        """
        Sort triangles.

        Parameters
        ----------
        fun : fun
            Key function.
        """

        self.items.sort(key=fun)

    #-----------------------------------------------------------------------------------------------

    def is_conflict(self, t):
        """
        Check is there conflict with triangle.

        Parameters
        ----------
        t : Triangle
            Triangle.

        Returns
        -------
        bool
            True - if there is conflict,
            False - otherwise.
        """

        for ti in self.items:
            if ti.is_conflict(t):
                return True

        return False

    #-----------------------------------------------------------------------------------------------

    def construct_mosaic(self):
        """
        Construct mosaic.
        """

        # First sort triangles from smallest.
        self.sort(fun=lambda t: t.area2())

        ts = self.items
        self.items = []

        # Try to add triangles without conflicts.
        for ti in ts:
            if not self.is_conflict(ti):
                self.add(ti)

#===================================================================================================

class Intersection:
    """
    Intersection of two geometrical objects.
    """

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def point_point(p1, p2):
        """
        Intersection of two points.

        Parameters
        ----------
        p1 : Point
            First point.
        p2 : Point.
            Secod point.

        Returns
        -------
        None
            No intersection.
        Point
            Equal points.
        """

        # Points may ne equal or not.
        if p1 == p2:
            return p1
        else:
            return None

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def point_line(p, ln):
        """
        Intersection of point and line.

        Parameters
        ----------
        p : Point
            Point.
        ln : Line
            Line.

        Returns
        -------
        None
            No intersection.
        Point
            Point is on line.
        """

        # Point is on line or there is no intersection.
        if ln.is_have_point(p):
            return p
        else:
            return None

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def point_segment(p, s):
        """
        Intersection of point and segment.

        Parameters
        ----------
        p : Point
            Point.
        s : Segment
            Segment.

        Returns
        -------
        None
            No intersection.
        Point
            Point is on segment.
        """

        # Point is on segment or there is no intersection.
        if s.is_have_point(p):
            return p
        else:
            return None

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def point_plane(p, pl):
        """
        Intersection of point and plane.

        Parameters
        ----------
        p : Point
            Point.
        pl : Plane
            Plane.

        Returns
        -------
        None
            No intersection.
        Point
            Point is on plane.
        """

        # Point is on plane or there is no intersection.
        if pl.is_have_point(p):
            return p
        else:
            return None

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def point_triangle(p, t):
        """
        Intersection of point and triangle.

        Parameters
        ----------
        p : Point
            Point.
        t : Triangle
            Triangle.

        Returns
        -------
        None
            No intersection.
        Point
            Point is inside triangle.
        """

        # Point is in triangle or there is no intersection.
        if t.is_have_point(p):
            return p
        else:
            return None

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def line_point(ln, p):
        """
        Intersection of line and point.

        Parameters
        ----------
        ln : Line
            Line.
        p : Point
            Point.

        Returns
        -------
        None
            No intersection.
        Point
            Point is on line.
        """

        return Intersection.point_line(p, ln)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def line_line(ln1, ln2):
        """
        Intersection of two lines.

        Parameters
        ----------
        ln1 : Line
            First line.
        ln2 : Line
            Second line.

        Returns
        -------
        None
            No intersection.
        Point
            Single points intersection.
        Line
            The same line.
        """

        # First check for equal.
        if ln1 == ln2:
            return ln1

        # Check for parallel lines.
        if Vector.vector_product(ln1.v, ln2.v).is_null():
            return None

        # Extract coefficients.
        x1, y1, z1 = ln1.P0.x, ln1.P0.y, ln1.P0.z
        m1, n1, p1 = ln1.v.x, ln1.v.y, ln1.v.z
        x2, y2, z2 = ln2.P0.x, ln2.P0.y, ln2.P0.z
        m2, n2, p2 = ln2.v.x, ln2.v.y, ln2.v.z

        # Linear equations system.
        # x1 + t1 * m1 = x2 + t2 * m2
        # y1 + t1 * n1 = y2 + t2 * n2
        # z1 + t1 * p1 = z2 + t2 * p2
        # Move all members to left.
        # t1 * m1 - t2 * m2 + (x1 - x2) = 0
        # t1 * n1 - t2 * n2 + (y1 - y2) = 0
        # t1 * p1 - t2 * p2 + (z1 - z2) = 0
        m2, n2, p2 = -m2, -n2, -p2
        dx, dy, dz = x1 - x2, y1 - y2, z1 - z2

        # System in simple form.
        # m1 * t1 + m2 * t2 + dx = 0 // (1)
        # n1 * t1 + n2 * t2 + dy = 0 // (2)
        # p1 * t1 + p2 * t2 + dz = 0 // (3)

        # Function for solving each system of 3.
        def slv2(m1, m2, dx, n1, n2, dy):
            # (1) * n1 - (2) * m1
            #   (m2 * n1 - n2 * m1) * t2 + (dx * n1 - dy * m1) = 0
            #   t2 = (dy * m1 - dx * n1) / q
            # (1) * n2 - (2) * m2
            #   (m1 * n2 - n1 * m2) * t1 + (dx * n2 - dy * m2) = 0
            #   t1 = (dx * n2 - dy * m2) / q
            q = m2 * n1 - n2 * m1
            if q == 0:
                return None
            else:
                return  (dx * n2 - dy * m2) / q, (dy * m1 - dx * n1) / q

        # Try to solve system of equations.
        r = slv2(m1, m2, dx, n1, n2, dy)
        if r is None:
            r = slv2(m1, m2, dx, p1, p2, dz)
            if r is None:
                r = slv2(n1, n2, dy, p1, p2, dz)
        t1, t2 = r

        # Try all equations.
        if (m1 * t1 + m2 * t2 + dx == 0) \
            and (n1 * t1 + n2 * t2 + dy == 0) \
            and (p1 * t1 + p2 * t2 + dz == 0):
            return Point(x1 + t1 * m1, y1 + t1 * n1, z1 + t1 * p1)
        else:
            return None

    #--------------------------------------------------------------------------------------

    @staticmethod
    def line_segment(ln, s):
        """
        Intersection of line and segment.

        Parameters
        ----------
        ln : Line
            Line.
        s : Segment
            Segment.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Segment is on line.
        """

        # Find intersection of two lines.
        r = Intersection.line_line(ln, s.line)

        # No intersection.
        if r is None:
            return None

        # Whole segment.
        if isinstance(r, Line):
            return s

        # Intersection of two lines is point.
        assert isinstance(r, Point)

        p = r

        # If point on segment then this point is intersection.
        if s.is_have_point(p):
            return p
        else:
            return None

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def line_plane(ln, pl):
        """
        Intersection of line and plane.

        Parameters
        ----------
        ln : Line
            Line.
        pl : Plane
            Plane.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Line
            Line is in plane.
        """

        x0, y0, z0, m, n, p = ln.P0.x, ln.P0.y, ln.P0.z, ln.v.x, ln.v.y, ln.v.z
        a, b, c, d = pl.a, pl.b, pl.c, pl.d

        # ax + by + cz + d = 0
        #   x = x0 + tm
        #   y = y0 + tn
        #   z = z0 + tp

        # One of m, n, p is not zero.

        if m != 0:

            # case 1. m != 0
            #   t = (x - x0) / m, y = y0 + (n / m)(x - x0), z = z0 + (p / m)(x - x0)
            #   ax + b (y0 + (n / m)(x - x0)) + c(z0 + (p / m)(x - x0)) + d = 0
            #   ax + b y0 + (bn / m)(x - x0) + c z0 + (cp / m)(x - x0) + d = 0
            #   ax + b y0 + (bn / m)x - (bn / m)x0 + c z0 + (cp / m)x - (cp / m)x0 + d = 0
            #   x (a + (bn + cp) / m) + b y0 + c z0 + d - ((bn + cp) / m) x0 = 0
            #   q = (bn + cp) / m
            #   x (a + q) + b y0 + c z0 + d - q x0 = 0
            #   x = (q x0 - b y0 - c z0 - d) / (a + q)
            q = (b * n + c * p) / m
            up, dn = q * x0 - b * y0 - c * z0 - d, a + q

            if dn == 0:
                if up == 0:
                    return ln
                else:
                    return None

            x = up / dn
            t = (x - x0) / m
            y = y0 + t * n
            z = z0 + t * p

            return Point(x, y, z)

        elif n != 0:

            # case 2. n != 0
            #   t = (y - y0) / n, x = x0 + (m / n)(y - y0), z = z0 + (p / n)(y - y0)
            #   a (x0 + (m / n)(y - y0)) + by + c (z0 + (p / n)(y - y0)) + d = 0
            #   a x0 + (am / n)(y - y0) + by + c z0 + (cp / n)(y - y0) + d = 0
            #   a x0 + (am / n)y - (am / n)y0 + by + c z0 + (cp / n)y - (cp / n)y0 + d = 0
            #   y ( b + (am + cp) / n) + a x0 + c z0 + d - ((am + cp) / n) y0 = 0
            #   q = (am + cp) / n
            #   y (b + q) + a x0 + c z0 + d - q y0 = 0
            #   y = (q y0 - a x0 - c z0 - d) / (b + q)
            q = (a * m + c * p) / n
            up, dn = q * y0 - a * x0 - c * z0 - d, b + q

            if dn == 0:
                if up == 0:
                    return ln
                else:
                    return None

            y = up / dn
            t = (y - y0) / n
            x = x0 + t * m
            z = z0 + t * p

            return Point(x, y, z)

        else:

            assert p != 0

            # case 3. p != 0
            #   t = (z - z0) / p, x = x0 + (m / p)(z - z0), y = y0 + (n / p)(z - z0)
            #   a (x0 + (m / p)(z - z0)) + b (y0 + (n / p)(z - z0)) + cz + d = 0
            #   a x0 + (am / p)(z - z0) + b y0 + (bn / p)(z - z0) + cz + d = 0
            #   a x0 + (am / p)z + (am / p)z0 + b y0 + (bn / p)z + (bn / p)z0 + cz + d = 0
            #   z (c + (am + bn) / p) + a x0 + b y 0 + d - ((am + bn) / p) z0 = 0
            #   q = (am + bn) / p
            #   z (c + q) + a x0 + b y0 + d - q z0 = 0
            #   z = (q z0 - a x0 - b y0 - d) / (c + q)
            q = (a * m + b * n) / p
            up, dn = q * z0 - a * x0 - b * y0 - d, c + q

            if dn == 0:
                if up == 0:
                    return ln
                else:
                    return None

            z = up / dn
            t = (z - z0) / p
            x = x0 + t * m
            y = y0 + t * n

            return Point(x, y, z)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def line_triangle(ln, t):
        """
        Intersection of line and triangle.

        Parameters
        ----------
        ln : Line
            Line.
        t : Triangle
            Triangle.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Intersection by segment.
        """

        raise Exception('not implemented')

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def segment_point(s, p):
        """
        Intersection of segment and point.

        Parameters
        ----------
        s : Segment
            Segment.
        p : Point
            Point.

        Returns
        -------
        None
            No intersection.
        Point
            Point is on segment.
        """

        return Intersection.point_segment(p, s)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def segment_line(s, ln):
        """
        Intersection of segment and line.

        Parameters
        ----------
        s : Segment
            Segment.
        ln : Line
            Line.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Segment is on line.
        """

        return Intersection.line_segment(ln, s)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def segment_segment(s1, s2):
        """
        Intersection of two segments.

        Parameters
        ----------
        s1 : Segment
            Segment.
        s2 : Segment
            Segment.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Intersection by segment.
        """

        # Find intersection of two containing lines.
        r = Intersection.line_line(s1.line, s2.line)

        # No intersection.
        if r is None:
            return None

        # One point intersection.
        # Both segments must have it.
        # We know that point is placed on both lines.
        # So we have to check point for between ends check.
        if isinstance(r, Point):
            p = r
            if s1.is_point_between_ends(p) and s2.is_point_between_ends(p):
                return p
            else:
                return None

        # So both segments lie in one line.
        assert isinstance(r, Line)

        # Get ends.
        A1, B1 = s1.A, s1.B
        A2, B2 = s2.A, s2.B
        assert A1 <= B1
        assert A2 <= B2

        # No intersection check.
        # A1       B1       A2       B2
        # *--------*........*--------*
        # A2       B2       A1       B1
        #*---------*........*--------*
        if (A2 > B1) or (A1 > B2):
            return None

        # Get intersection segment.
        # A1       A2       B1       B2
        # *--------*========*--------*
        # A2       A1       B2       B1
        # *--------*========*--------*
        A, B = max(A1, A2), min(B1, B2)

        # Intersection can be segment or single point.
        if A == B:
            return A
        else:
            return Segment(A, B)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def segment_plane(s, pl):
        """
        Intersection of segment and plane.

        Parameters
        ----------
        s : Segment
            Segment.
        pl : Plane
            Plane.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        """

        # Check for no intersection.
        if pl.is_two_points_strong_on_one_side(s.A, s.B):
            return None

        # Check if both ends of segment lie in plane.
        if pl.is_have_point(s.A) and pl.is_have_point(s.B):
            return s

        # Intersection point is inside of segment.
        return Intersection.line_plane(s.line, pl)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def segment_triangle(s, t):
        """
        Intersection of segment and triangle.

        Parameters
        ----------
        s : Segment
            Segment.
        t : Triangle
            Triangle.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Intersection by segment.
        """

        # Find intersection of plane and segment.
        r = Intersection.segment_plane(s, t.plane)

        # No intersection segment with plane.
        if r is None:
            return None

        # Single point on plane.
        if isinstance(r, Point):
            p = r
            if t.is_have_point(p):
                return p
            else:
                return None

        # Segment is in triangle plane.
        assert isinstance(r, Segment)

        # Find intersections of all triangle sides with segment.
        r1 = Intersection.segment_segment(t.AB, s)
        r2 = Intersection.segment_segment(t.BC, s)
        r3 = Intersection.segment_segment(t.AC, s)

        # If one of intersections r1, r2, r3 is segment,
        # then this segment is result of triangle and segment intersection
        # and its whole triangle side.
        if isinstance(r1, Segment):
            return r1
        if isinstance(r2, Segment):
            return r2
        if isinstance(r3, Segment):
            return r3

        # Then we have to collect all points.
        ps = Points()
        if isinstance(r1, Point):
            ps.add_unique(r1)
        if isinstance(r2, Point):
            ps.add_unique(r2)
        if isinstance(r3, Point):
            ps.add_unique(r3)

        # Now result is None, Point or Segment.
        cnt = ps.count()

        # Both intersection points are on sides.
        # They form result segment.
        if cnt == 2:
            return Segment(ps[0], ps[1])

        # Only one intersection point on side.
        if cnt == 1:
            if t.is_have_point(s.A) and t.is_have_point(s.B):
                return s
            else:
                return ps[0]

        # No intersection points with triangle sides.
        assert cnt == 0

        # Result if whole segment or empty.
        if t.is_have_point(s.A) and t.is_have_point(s.B):
            return s
        else:
            return None

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def plane_point(pl, p):
        """
        Intersection of plane and segment.

        Parameters
        ----------
        pl : Plane
            Plane.
        p : Point
            Point.

        Returns
        -------
        None
            No intersection.
        Point
            Point is on plane.
        """

        return Intersection.point_plane(p, pl)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def plane_line(pl, ln):
        """
        Intersection of plane with line.

        Parameters
        ----------
        pl : Plane
            Plane.
        ln : Line
            Line.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Line
            Line is in plane.
        """

        return Intersection.line_plane(ln, pl)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def plane_segment(pl, s):
        """
        Intersection of plane and segment.

        Parameters
        ----------
        pl : Plane
            Plane.
        s : Segment
            Segment.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Segment is in plane.
        """

        return Intersection.segment_plane(s, pl)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def plane_plane(pl1, pl2):
        """
        Intersection of two planes.

        Parameters
        ----------
        pl1 : Plane
            First plane.
        pl2 : Plane
            Second plane.

        Returns
        -------
        None
            Parallel planes.
        Line
            Intersect planes.
        """

        raise Exception('not implemented')

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def plane_triangle(pl, t):
        """
        Intersection of plane and triangle.

        Parameters
        ----------
        pl : Plane
            Plane.
        t : Triangle
            Triangle.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Segment intersection.
        Triangle
            Triangle is in plane.
        """

        raise Exception('not implemented')

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def triangle_point(t, p):
        """
        Intersection of triangle with point.

        Parameters
        ----------
        t : Triangle
            Triangle.
        p : Point
            Point.

        Returns
        -------
        None
            No intersection.
        Point
            Point is in triangle.
        """

        return Intersection.point_triangle(p, t)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def triangle_line(t, ln):
        """
        Intersection of triangle and line.

        Parameters
        ----------
        t : Triangle
            Triangle.
        ln : Line
            Line.

        Returns
        -------
        None
            No intersection.
        Point
            Single point intersection.
        Segment
            Intersection by segment.
        """

        return Intersection.line_triangle(ln, t)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def triangle_segment(t, s):
        """
        Intersection of triangle with segment.

        Parameters
        ----------
        t : Triangle
            Triangle.
        s : Segment
            Segment.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Intersection by segment.
        """

        return Intersection.segment_triangle(s, t)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def triangle_plane(t, pl):
        """
        Intersection triangle with plane.

        Parameters
        ----------
        t : Triangle
            Triangle.
        pl : Plane
            Plane.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Segment intersection.
        Triangle
            Triangle is in plane.
        """

        return Intersection.plane_triangle(pl, t)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def triangle_triangle(t1, t2):
        """
        Intersection of two triangles.

        Parameters
        ----------
        t1 : Triangle
            First triangle.
        t2 : Triangle
            Second triangle.

        Returns
        -------
        None
            No intersection.
        Point
            Single point of intersection.
        Segment
            Intersection by segment.
        Triangle
            Triangle intersection.
        Points
            Convex figure of intersection
            4 - convex quadrangle,
            5 - convex pentagon,
            6 - convex hexagon.
        """

        # Fast check for None intersection.
        if t1.plane.is_three_points_strong_on_one_side(t2.A, t2.B, t2.C):
            return None
        if t2.plane.is_three_points_strong_on_one_side(t1.A, t1.B, t1.C):
            return None

        ps = Points()

        for s in t1.sides:
            r = Intersection.segment_triangle(s, t2)
            if isinstance(r, Point):
                ps.add_unique(r)
            elif isinstance(r, Segment):
                ps.add_unique(r.A)
                ps.add_unique(r.B)

        for s in t2.sides:
            r = Intersection.segment_triangle(s, t1)
            if isinstance(r, Point):
                ps.add_unique(r)
            elif isinstance(r, Segment):
                ps.add_unique(r.A)
                ps.add_unique(r.B)

        cnt = ps.count()

        # Check all cases.
        if cnt == 0:
            return None
        elif cnt == 1:
            return ps[0]
        elif cnt == 2:
            return Segment(ps[0], ps[1])
        elif cnt == 3:
            ss = Segments()
            ss.adds([Segment(ps[0], ps[1]), Segment(ps[1], ps[2]), Segment(ps[0], ps[1])])
            return ss
        else:
            raise Exception('Intersection.triangle_triangle NOT IMPLEMENTED')

#===================================================================================================
# Global objects.
#===================================================================================================

# Base points.
O = Point(Fr(0), Fr(0), Fr(0))
X = Point(Fr(1), Fr(0), Fr(0))
Y = Point(Fr(0), Fr(1), Fr(0))
Z = Point(Fr(0), Fr(0), Fr(1))

# Base lines.
OX = Line.from_points(O, X)
OY = Line.from_points(O, Y)
OZ = Line.from_points(O, Z)

# Base planes.
OXY = Plane.from_points(O, X, Y)
OYZ = Plane.from_points(O, Y, Z)
OXZ = Plane.from_points(O, X, Z)

#===================================================================================================
# Global functions.
#===================================================================================================

def triangles_list_box(ts, eps):
    """
    Box of triangles list.

    Parameters
    ----------
    ts : [Triangle]
        List of triangles.
    eps : float
        Epsilon.

    Returns
    -------
    geom3d.Box
        Box.
    """

    lo, hi = float('inf'), float('-inf')
    xlo, xhi, ylo, yhi, zlo, zhi = lo, hi, lo, hi, lo, hi

    # Check all triangles.
    for t in ts:
        ps = t.points
        xlo, xhi = min(xlo, min([p.x for p in ps])), max(xhi, max([p.x for p in ps]))
        ylo, yhi = min(ylo, min([p.y for p in ps])), max(yhi, max([p.y for p in ps]))
        zlo, zhi = min(zlo, min([p.z for p in ps])), max(zhi, max([p.z for p in ps]))

    return geom3d.Box(geom1d.Segment(float(xlo), float(xhi)),
                      geom1d.Segment(float(ylo), float(yhi)),
                      geom1d.Segment(float(zlo), float(zhi)), eps)

#---------------------------------------------------------------------------------------------------

def split_triangles_list(ts, d, v):
    """
    Split triangles list into two.

    Parameters
    ----------
    ts : [(Triangle, int)]
        List of pairs - triangle / index.
    d : str
        Direction ('x', 'y', 'z').
    v : float
        Split value.

    Returns
    -------
    [(Triangle, int)], [(Triangle, int)]
        Splitted structure.
    """

    ts_lo, ts_hi = [], []

    def anyloeq(v1, v2, v3, v):
        return (v1 <= v) or (v2 <= v) or (v3 <= v)

    def anyhieq(v1, v2, v3, v):
        return (v1 >= v) or (v2 >= v) or (v3 >= v)

    if d == 'x':
        for (t, i) in ts:
            v1, v2, v3 = t.A.x, t.B.x, t.C.x
            if anyloeq(v1, v2, v3, v):
                ts_lo.append((t, i))
            if anyhieq(v1, v2, v3, v):
                ts_hi.append((t, i))
    elif d == 'y':
        for (t, i) in ts:
            v1, v2, v3 = t.A.y, t.B.y, t.C.y
            if anyloeq(v1, v2, v3, v):
                ts_lo.append((t, i))
            if anyhieq(v1, v2, v3, v):
                ts_hi.append((t, i))
    else:
        if d != 'z':
            raise Exception('geom3d_rat:split_triangles_Liist: wrong split direction.')
        for (t, i) in ts:
            v1, v2, v3 = t.A.z, t.B.z, t.C.z
            if anyloeq(v1, v2, v3, v):
                ts_lo.append((t, i))
            if anyhieq(v1, v2, v3, v):
                ts_hi.append((t, i))

    return ts_lo, ts_hi

#---------------------------------------------------------------------------------------------------

def best_split_triangles_list(ts, vx, vy, vz):
    """
    Check best split.

    Parameters
    ----------
    ts : [(Triangle, int)]
        List of pairs - triangle / index.
    vx : float
        Value for X split.
    vy : float
        Value for Y split.
    vz : float
        Value for Z split.

    Returns
    -------
    (str, (Triangle, int), (Triangle, int)) | None
        Splitted structure with direction,
        or None - if split is not possible.
    """

    ts_xlo, ts_xhi = split_triangles_list(ts, 'x', vx)
    ts_ylo, ts_yhi = split_triangles_list(ts, 'y', vy)
    ts_zlo, ts_zhi = split_triangles_list(ts, 'z', vz)
    n = len(ts)
    ts_xlon, ts_xhin = len(ts_xlo), len(ts_xhi)
    ts_ylon, ts_yhin = len(ts_ylo), len(ts_yhi)
    ts_zlon, ts_zhin = len(ts_zlo), len(ts_zhi)
    dx, dy, dz = ts_xlon + ts_xhin - n, ts_ylon + ts_yhin - n, ts_zlon + ts_zhin - n

    def is_valid_sec(lo, hi, n):
        return (lo != 0) and (lo != n) and (hi != 0) and (hi != n)

    if (dx <= dy) and (dx <= dz) and is_valid_sec(ts_xlon, ts_xhin, n):
        # x - is the best
        return 'x', ts_xlo, ts_xhi
    elif (dy <= dz) and is_valid_sec(ts_ylon, ts_yhin, n):
        # y - is the best
        return 'y', ts_ylo, ts_yhi
    elif is_valid_sec(ts_zlon, ts_zhin, n):
        return 'z', ts_zlo, ts_zhi
    else:
        return None

#---------------------------------------------------------------------------------------------------

def split_bvh_tree_with_rat_triangles(bvh):
    """
    Split BVH tree with rat triangles.

    Parameters
    ----------
    bvh : bvh_tree.BVHTRee
        Tree.
    """

    if bvh.children:
        raise Exception('geom3d_rat.split_bvh_tree_with_rat_triangles: '
                        'can not split BVH tree with children.')

    [vx, vy, vz] = bvh.box.mids()
    bs = best_split_triangles_list(bvh.data, vx, vy, vz)

    if bs is None:
        return

    d, lo, hi = bs
    bvh.split(d)
    bvh.children[0].data = lo
    bvh.children[1].data = hi

    # Recursive split.
    split_bvh_tree_with_rat_triangles(bvh.children[0])
    split_bvh_tree_with_rat_triangles(bvh.children[1])

#---------------------------------------------------------------------------------------------------

def rat_triangles_bvh_tree(ts):
    """
    Create rat triangles BVH tree.

    Parameters
    ----------
    ts : [Triangle()]
        List of triangles.

    Returns
    -------
    bvh_tree.BVHTRee
        BVH tree.
    """

    # Wrap triangles.
    tss = [(ts[i], i) for i in range(len(ts))]

    # Create root of BVH.
    bvh = bvh_tree.BVHTree(triangles_list_box(ts, 0.0))
    bvh.data = tss
    split_bvh_tree_with_rat_triangles(bvh)

    return bvh

#---------------------------------------------------------------------------------------------------

def print_statistics():
    """
    Print statistics.
    """

    print(f'geom3d_rat created objects: {Point.counter} points, '
          f'{Line.counter} lines, {Segment.counter} segments, '
          f'{Plane.counter} planes, {Triangle.counter} triangles')

#===================================================================================================

def test():
    """
    Tests.
    """

    # objects.
    XY = Line.from_points(X, Y)
    XYZ = Plane.from_points(X, Y, Z)

    # Check compare points.
    assert Point(Fr(1), Fr(0), Fr(0)) < Point(Fr(2), Fr(0), Fr(0))
    assert Point(Fr(1), Fr(2), Fr(0)) > Point(Fr(1), Fr(1), Fr(0))
    assert Point(Fr(0), Fr(0), Fr(4)) < Point(Fr(0), Fr(0), Fr(6))
    assert Point(Fr(2), Fr(3), Fr(4)) == Point(Fr(2), Fr(3), Fr(4))

    # Check line have points.
    assert Line(Fr(0), Fr(0), Fr(0), Fr(1), Fr(1), Fr(1)).is_have_point(Point(Fr(2), Fr(2), Fr(2)))

    # Perpendicular planes.
    assert OXY.is_perpendicular_with_plane(OYZ)
    assert OXY.is_perpendicular_with_plane(OXZ)
    assert OYZ.is_perpendicular_with_plane(OXZ)

    # Intersect plane with line.
    assert Intersection.line_plane(OX, XYZ) == X
    assert Intersection.line_plane(OY, XYZ) == Y
    assert Intersection.line_plane(OZ, XYZ) == Z
    assert Intersection.line_plane(OX, OXY) == OX
    assert Intersection.line_plane(OY, OYZ) == OY
    assert Intersection.line_plane(OZ, OXZ) == OZ

    # Check equal lines.
    assert OX == Line.from_points(O, Point(Fr(2), Fr(0), Fr(0)))

    # Check point on segment.
    SOX = Segment(O, X)
    assert SOX.is_have_point(O)
    assert SOX.is_have_point(Point(Fr(1, 2), Fr(0), Fr(0)))
    assert not SOX.is_have_point(Y)
    assert not SOX.is_have_point(Point(Fr(3, 2), Fr(0), Fr(0)))

    # Check point on triangle.
    TXYZ = Triangle(X, Y, Z)
    assert TXYZ.is_have_point(Point(Fr(1, 3), Fr(1, 3), Fr(1, 3)))
    assert TXYZ.is_have_point(Point(Fr(2, 5), Fr(2, 5), Fr(1, 5)))
    assert not TXYZ.is_have_point(Point(Fr(2, 3), Fr(2, 3), Fr(-1, 3)))
    assert not TXYZ.is_have_point(Point(Fr(-1), Fr(-1), Fr(-1)))

    # Check lines intersections.
    assert Intersection.line_line(OX, OY) == O
    assert Intersection.line_line(OX, OZ) == O

    # Check segments intersections.
    p1 = Point(Fr(-18789, 1000000), Fr(13811, 1000000), Fr(1, 100))
    p2 = Point(Fr(-9207, 500000), Fr(1841, 125000), Fr(1, 100))
    p3 = Point(Fr(-13092327, 705500000), Fr(5071507, 352750000), Fr(1, 100))
    s1 = Segment(p1, p2)
    s2 = Segment(p2, p3)
    assert s1.is_conflict(s2)

    # Segment and triangle intersection.
    s = Segment(Point(Fr(-1), Fr(0), Fr(0)), Point(Fr(0), Fr(1), Fr(0)))
    t = Triangle(Point(Fr(-1), Fr(0), Fr(0)),
                  Point(Fr(1), Fr(0), Fr(0)),
                  Point(Fr(0), Fr(2), Fr(0)))
    r = Intersection.segment_triangle(s, t)
    assert r == s
    s = Segment(Point(Fr(0), Fr(1, 2), Fr(0)), Point(Fr(0), Fr(1), Fr(0)))
    t = Triangle(Point(Fr(-1), Fr(0), Fr(0)),
                  Point(Fr(1), Fr(0), Fr(0)),
                  Point(Fr(0), Fr(2), Fr(0)))
    r = Intersection.segment_triangle(s, t)
    assert r == s

    # Intersection of two triangles.
    t1 = Triangle(Point(Fr(0), Fr(1), Fr(0)),
                  Point(Fr(1), Fr(-1), Fr(0)),
                  Point(Fr(0), Fr(-1), Fr(1)))
    t2 = Triangle(Point(Fr(0), Fr(0), Fr(1, 10)),
                  Point(Fr(1), Fr(2), Fr(1, 10)),
                  Point(Fr(0), Fr(2), Fr(11, 10)))
    r = Intersection.triangle_triangle(t1, t2)
    assert r == Segment(Point(Fr(0), Fr(2, 5), Fr(3, 10)), Point(Fr(1, 5), Fr(2, 5), Fr(1, 10)))
    t1 = Triangle(Point(Fr(-1), Fr(0), Fr(0)),
                  Point(Fr(1), Fr(0), Fr(0)),
                  Point(Fr(0), Fr(1), Fr(0)))
    t2 = Triangle(Point(Fr(-1), Fr(0), Fr(0)),
                  Point(Fr(1), Fr(0), Fr(0)),
                  Point(Fr(0), Fr(2), Fr(0)))
    r = Intersection.triangle_triangle(t1, t2)
    assert isinstance(r, Segments)
    assert r.count() == 3

    # Segment split.
    s = Segment(Point(Fr(0), Fr(0), Fr(0)), Point(Fr(1), Fr(0), Fr(0)))
    ps = Points()
    ps.add_unique(Point(Fr(-1), Fr(0), Fr(0)))
    assert s.split(ps) == s
    ps.add_unique(Point(Fr(0), Fr(0), Fr(0)))
    assert s.split(ps) == s
    ps.add_unique(Point(Fr(1, 2), Fr(0), Fr(0)))
    r = s.split(ps)
    assert (len(r) == 2)
    assert r[0] == Segment(Point(Fr(0), Fr(0), Fr(0)), Point(Fr(1, 2), Fr(0), Fr(0)))
    assert r[1] == Segment(Point(Fr(1, 2), Fr(0), Fr(0)), Point(Fr(1), Fr(0), Fr(0)))
    ps.add_unique(Point(Fr(3, 4), Fr(0), Fr(0)))
    r = s.split(ps)
    assert (len(r) == 3)
    assert r[0] == Segment(Point(Fr(0), Fr(0), Fr(0)), Point(Fr(1, 2), Fr(0), Fr(0)))
    assert r[1] == Segment(Point(Fr(1, 2), Fr(0), Fr(0)), Point(Fr(3, 4), Fr(0), Fr(0)))
    assert r[2] == Segment(Point(Fr(3, 4), Fr(0), Fr(0)), Point(Fr(1), Fr(0), Fr(0)))
    ps.add_unique(Point(Fr(2), Fr(0), Fr(0)))
    r = s.split(ps)
    assert (len(r) == 3)

#---------------------------------------------------------------------------------------------------

def test_triangulation(N):
    """
    Test triangulation.

    Parameters
    ----------
    N : int
        Test number.
    """

    if N == 1:
        # One vertical segment inside.
        t = Triangle(Point(Fr(1), Fr(0), Fr(0)),
                     Point(Fr(-1), Fr(0), Fr(0)),
                     Point(Fr(0), Fr(2), Fr(0)))
        pas = PointsAndSegments()
        pas.add_unique_segment(Segment(Point(Fr(0), Fr(1, 2), Fr(0)),
                                       Point(Fr(0), Fr(1), Fr(0))))
    elif N == 2:
        # One side segment.
        t = Triangle(Point(Fr(1), Fr(0), Fr(0)),
                     Point(Fr(-1), Fr(0), Fr(0)),
                     Point(Fr(0), Fr(2), Fr(0)))
        pas = PointsAndSegments()
        pas.add_unique_segment(Segment(Point(Fr(-1, 4), Fr(1, 2), Fr(0)),
                                       Point(Fr(1, 4), Fr(1), Fr(0))))
    elif N == 3:
        # Two not conflicted segments.
        t = Triangle(Point(Fr(1), Fr(0), Fr(0)),
                     Point(Fr(-1), Fr(0), Fr(0)),
                     Point(Fr(0), Fr(2), Fr(0)))
        pas = PointsAndSegments()
        pas.add_unique_segment(Segment(Point(Fr(-1, 2), Fr(1, 10), Fr(0)),
                                       Point(Fr(0), Fr(1), Fr(0))))
        pas.add_unique_segment(Segment(Point(Fr(1, 2), Fr(3, 10), Fr(0)),
                                       Point(Fr(1, 2), Fr(1, 10), Fr(0))))
    elif N == 4:
        # Whole triangle inside.
        t = Triangle(Point(Fr(1), Fr(0), Fr(0)),
                     Point(Fr(-1), Fr(0), Fr(0)),
                     Point(Fr(0), Fr(2), Fr(0)))
        pas = PointsAndSegments()
        pas.add_unique_segment(Segment(Point(Fr(-1, 2), Fr(1, 10), Fr(0)),
                                       Point(Fr(1, 2), Fr(1, 10), Fr(0))))
        pas.add_unique_segment(Segment(Point(Fr(-1, 2), Fr(1, 10), Fr(0)),
                                       Point(Fr(0), Fr(1), Fr(0))))
        pas.add_unique_segment(Segment(Point(Fr(1, 2), Fr(1, 10), Fr(0)),
                                       Point(Fr(0), Fr(1), Fr(0))))
    elif N == 5:
        # One vertical segment with one end on bottom side of triangle.
        t = Triangle(Point(Fr(1), Fr(0), Fr(0)),
                     Point(Fr(-1), Fr(0), Fr(0)),
                     Point(Fr(0), Fr(2), Fr(0)))
        pas = PointsAndSegments()
        pas.add_unique_segment(Segment(Point(Fr(0), Fr(0), Fr(0)),
                                       Point(Fr(0), Fr(1), Fr(0))))
    elif N == 6:
        # Segment lies on two sides.
        t = Triangle(Point(Fr(1), Fr(0), Fr(0)),
                     Point(Fr(-1), Fr(0), Fr(0)),
                     Point(Fr(0), Fr(2), Fr(0)))
        pas = PointsAndSegments()
        pas.add_unique_segment(Segment(Point(Fr(1, 2), Fr(1), Fr(0)),
                                       Point(Fr(-1, 2), Fr(1), Fr(0))))
    elif N == 7:
        # issue #5
        # https://github.com/r-aax/rumatha/issues/5
        # Intersection from cyl_1_right pseudo 3D profile.
        t = Triangle(Point(Fr(-9381, 500000), Fr(14131, 1000000), Fr(0)),
                     Point(Fr(-18789, 1000000), Fr(13811, 1000000), Fr(1, 100)),
                     Point(Fr(-9207, 500000), Fr(1841, 125000), Fr(1, 100)))
        pas = PointsAndSegments()
        p1 = Point(Fr(-84787749, 4553500000), Fr(4090753, 284593750), Fr(3707, 910700))
        p2 = Point(Fr(-13092327, 705500000), Fr(5071507, 352750000), Fr(1, 100))
        assert t.is_have_point_on_sides(p1) and t.is_have_point_on_sides(p2)
        s = Segment(p1, p2)
        pas.add_unique_segment(s)
    elif N == 8:
        # Many trajectories.
        t = Triangle(Point(Fr(1), Fr(1), Fr(0)),
                     Point(Fr(11), Fr(1), Fr(0)),
                     Point(Fr(6), Fr(8), Fr(0)))
        pas = PointsAndSegments()
        # 1st trajectory
        pas.add_unique_segment(Segment(Point(Fr(2), Fr(12, 5), Fr(0)),
                                       Point(Fr(6), Fr(2), Fr(0))))
        pas.add_unique_segment(Segment(Point(Fr(6), Fr(2), Fr(0)),
                                       Point(Fr(8), Fr(26, 5), Fr(0))))
        # 2nd trajectory
        pas.add_unique_segment(Segment(Point(Fr(8), Fr(1), Fr(0)),
                                       Point(Fr(19, 2), Fr(31, 10), Fr(0))))
        # 3rd trajectory
        pas.add_unique_segment(Segment(Point(Fr(7, 2), Fr(1), Fr(0)),
                                       Point(Fr(3), Fr(3), Fr(0))))
        pas.add_unique_segment(Segment(Point(Fr(3), Fr(3), Fr(0)),
                                       Point(Fr(5), Fr(11, 2), Fr(0))))
        pas.add_unique_segment(Segment(Point(Fr(5), Fr(11, 2), Fr(0)),
                                       Point(Fr(7), Fr(33, 5), Fr(0))))
        # 4th trajectory
        pas.add_unique_segment(Segment(Point(Fr(4), Fr(26, 5), Fr(0)),
                                       Point(Fr(5), Fr(7, 2), Fr(0))))
        pas.add_unique_segment(Segment(Point(Fr(5), Fr(7, 2), Fr(0)),
                                       Point(Fr(6), Fr(5), Fr(0))))
        pas.add_unique_segment(Segment(Point(Fr(6), Fr(5), Fr(0)),
                                       Point(Fr(8), Fr(3), Fr(0))))
        pas.add_unique_segment(Segment(Point(Fr(8), Fr(3), Fr(0)),
                                       Point(Fr(7), Fr(1), Fr(0))))
    elif N == 9:
        # With points
        t = Triangle(Point(Fr(1), Fr(0), Fr(0)),
                     Point(Fr(-1), Fr(0), Fr(0)),
                     Point(Fr(0), Fr(2), Fr(0)))
        pas = PointsAndSegments()
        pas.add_unique_point(Point(Fr(0), Fr(5, 4), Fr(0)))
        pas.add_unique_point(Point(Fr(1, 3), Fr(1, 3), Fr(0)))
        pas.add_unique_point(Point(Fr(-1, 3), Fr(1, 3), Fr(0)))
    elif N == 10:
        # With points and segment.
        t = Triangle(Point(Fr(1), Fr(0), Fr(0)),
                     Point(Fr(-1), Fr(0), Fr(0)),
                     Point(Fr(0), Fr(2), Fr(0)))
        pas = PointsAndSegments()
        pas.add_unique_segment(Segment(Point(Fr(1, 2), Fr(3, 5), Fr(0)),
                                       Point(Fr(-1, 2), Fr(3, 5), Fr(0))))
        pas.add_unique_point(Point(Fr(1, 2), Fr(4, 5), Fr(0)))
        pas.add_unique_point(Point(Fr(1, 2), Fr(2, 5), Fr(0)))
        pas.add_unique_point(Point(Fr(-1, 2), Fr(4, 5), Fr(0)))
        pas.add_unique_point(Point(Fr(-1, 2), Fr(2, 5), Fr(0)))
        pas.add_unique_point(Point(Fr(0), Fr(4, 5), Fr(0)))
        pas.add_unique_point(Point(Fr(0), Fr(2, 5), Fr(0)))
    else:
        assert False

    #
    # Final processing.
    #

    # Add triangle to pas.
    pas.adds_unique_segments(t.sides)
    ts = None

    # Triangulate.
    ts = t.triangulate(pas)

    # Draw canvas.
    pas.draw(plt, color='orange', linewidth='5', size=60)

    if not ts is None:
        for ti in ts:
            ti.draw(plt, color='black', linewidth='1', size=0)
    plt.show()

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    start = time.time()
    test()
    #test_triangulation(N=10)
    print_statistics()
    print(f'total time : {time.time() - start}')

#===================================================================================================
