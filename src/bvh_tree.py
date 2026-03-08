"""
BVH-tree for geometrical objects processing speedup.
"""

import geom3d

#===================================================================================================

class BVHTree:
    """
    BVH tree.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, box):
        """
        Create BVH-tree.

        Parameters
        ----------
        box : Box
            Box.
        """

        self.box = box
        self.children = []

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def from_floats(xlo, xhi, ylo, yhi, zlo, zhi, eps):
        """

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
        eps : float
            Epsilon.

        Returns
        -------
        BVHTree
            Created tree.
        """

        # Convert bounds to floats.
        xlo, xhi, ylo, yhi, zlo, zhi = float(xlo), float(xhi), \
                                       float(ylo), float(yhi), float(zlo), float(zhi)

        box = geom3d.Box.from_floats(xlo, xhi, ylo, yhi, zlo, zhi, eps)

        return BVHTree(box)

    #-----------------------------------------------------------------------------------------------

    def print(self, level=0):
        """
        Print BVH-tree information.

        Parameters
        ----------
        level : int
            Level.
        """

        off = ' ' * level

        print(f'{off} l{level}: box {self.box}')

        for ch in self.children:
            ch.print(level + 1)

    #-----------------------------------------------------------------------------------------------

    def split(self, d):
        """
        Split BVH tree by direction.

        Parameters
        ----------
        d : str
            Direction ('x', 'y', 'z').
        """

        if self.children:
            raise Exception('bvh_tree:BVHTree.split: can not split tree which has children.')

        b1, b2 = self.box.split(d)
        self.children.append(BVHTree(b1))
        self.children.append(BVHTree(b2))

#===================================================================================================

if __name__ == '__main__':
    b = BVHTree.from_floats(0, 100, 0, 100, 0, 100, 0.001)
    b.split('x')
    b.children[0].split('y')
    b.children[0].children[1].split('z')
    b.print()

#===================================================================================================
