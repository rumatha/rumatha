"""
Cutting blocks.
"""

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import random

matplotlib.rcParams.update({'font.size': 22})

#===================================================================================================

# Find minimal cuts to extract from block n x m x k (n >= m >= k)
# part of size t.

MIN_CUTS_FOR_EXTRACT_PART_ON = False

if MIN_CUTS_FOR_EXTRACT_PART_ON:
    # Global mem (we can process blocks with n <= 10, m <= 10, k <= 10).
    min_cuts_for_extract_part_n = 11
    min_cuts_for_extract_part_shape = (min_cuts_for_extract_part_n,
                                       min_cuts_for_extract_part_n,
                                       min_cuts_for_extract_part_n,
                                       min_cuts_for_extract_part_n**3)
    min_cuts_for_extract_part_mem = 0 - np.ones(min_cuts_for_extract_part_shape, dtype=int)

#---------------------------------------------------------------------------------------------------

def min_cuts_for_extract_part_1d(n, t):
    """
    Find minimum cuts count to extract from block n
    part of size t.

    Parameters
    ----------
    n : int
        Block side size.
    t : int
        Target part size.

    Returns
    -------
    int
        Minimum cuts count.
    """

    if t > n:
        return math.inf

    # Use memory.
    if min_cuts_for_extract_part_mem[n][1][1][t] >= 0:
        return min_cuts_for_extract_part_mem[n][1][1][t]

    if (t == 0) or (t == n):
        r = 0
    else:
        r = 1

    # Use memory.
    min_cuts_for_extract_part_mem[n][1][1][t] = r

    return r

#---------------------------------------------------------------------------------------------------

def min_cuts_for_extract_part_2d(n, m, t):
    """
    Find minimal cut count to extract from block n x m
    part of size t.

    Parameters
    ----------
    n : int
        Block side size.
    m : int
        Block side size.
    t : int
        Target part size.

    Returns
    -------
    int
        Minimal cuts count.
    """

    s = n * m

    if t > s:
        return math.inf

    # Normalize n >= m, t <= s / 2.
    if not (n >= m):
        n, m = m, n
    if t > s // 2:
        t = s - t

    # Use memory.
    if min_cuts_for_extract_part_mem[n][m][1][t] >= 0:
        return min_cuts_for_extract_part_mem[n][m][1][t]

    if (t == 0) or (t == s):
        r = 0
    elif n == 1:
        r = min_cuts_for_extract_part_1d(m, t)
    elif m == 1:
        r = min_cuts_for_extract_part_1d(n, t)
    else:
        r = math.inf
        for t1 in range(t):
            for n1 in range(1, n):
                r = min(r,
                        min_cuts_for_extract_part_2d(n1, m, t1) \
                        + min_cuts_for_extract_part_2d(n - n1, m, t - t1))
            for m1 in range(1, m):
                r = min(r,
                        min_cuts_for_extract_part_2d(n, m1, t1) \
                        + min_cuts_for_extract_part_2d(n, m - m1, t - t1))
        r = r + 1

    # Use memory.
    min_cuts_for_extract_part_mem[n][m][1][t] = r

    return r

#---------------------------------------------------------------------------------------------------

def min_cuts_for_extract_part_3d(n, m, k, t):
    """
    Find minimal cuts count to extract from block n x m x k
    part of size t.

    Parameters
    ----------
    n : int
        Block side size.
    m : int
        Block side size.
    k : int
        Block side size.
    t : int
        Target part size.

    Returns
    -------
    int
        Minimal cuts count.
    """

    s = n * m * k

    if t > s:
        return math.inf

    # Normalize n >= m >= k, t <= s / 2
    if not ((n >= m) and (m >= k)):
        a = [n, m, k]
        a.sort()
        [k, m, n] = a
    if t > s // 2:
        t = s - t

    # Use memory.
    if min_cuts_for_extract_part_mem[n][m][k][t] >= 0:
        return min_cuts_for_extract_part_mem[n][m][k][t]

    if (t == 0) or (t == s):
        r = 0
    elif n == 1:
        r = min_cuts_for_extract_part_2d(m, k, t)
    elif m == 1:
        r = min_cuts_for_extract_part_2d(n, k, t)
    elif k == 1:
        r = min_cuts_for_extract_part_2d(n, m, t)
    else:
        r = math.inf
        for t1 in range(t):
            for n1 in range(1, n):
                r = min(r,
                        min_cuts_for_extract_part_3d(n1, m, k, t1) \
                        + min_cuts_for_extract_part_3d(n - n1, m, k, t - t1))
            for m1 in range(1, m):
                r = min(r,
                        min_cuts_for_extract_part_3d(n, m1, k, t1) \
                        + min_cuts_for_extract_part_3d(n, m - m1, k, t - t1))
            for k1 in range(1, k):
                r = min(r,
                        min_cuts_for_extract_part_3d(n, m, k1, t1) \
                        + min_cuts_for_extract_part_3d(n, m, k - k1, t - t1))
        r = r + 1

    # Use memory.
    min_cuts_for_extract_part_mem[n][m][k][t] = r

    return r

#===================================================================================================

class Block:
    """
    Block of sizes a * b * c.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, a, b, c):
        """
        Block constructor.

        Parameters
        ----------
        a : int
            Side size.
        b : int
            Side size.
        c : int
            Side size.
        """

        # We need first dimension maximum.
        li = [a, b, c]
        li.sort(reverse=True)
        [a, b, c] = li

        self.a = a
        self.b = b
        self.c = c

    #-----------------------------------------------------------------------------------------------

    def copy(self):
        """
        Copy constructor.

        Returns
        -------
        Block
            Block copy.
        """

        return Block(self.a, self.b, self.c)

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def random(alo, ahi, blo, bhi, clo, chi):
        """
        Create random block.

        Parameters
        ----------
        alo : int
            Lo value of side a.
        ahi : int
            Hi value of side a.
        blo : int
            Lo value of side b.
        bhi : int
            Hi value of side b.
        clo : int
            Lo value of side c.
        chi : int
            Hi value of side c.

        Returns
        -------
        Block
            Random block.
        """

        return Block(random.randint(alo, ahi),
                     random.randint(blo, bhi),
                     random.randint(clo, chi))

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'B({self.a}-{self.b}-{self.c},{self.cells_count}c)'

    #-----------------------------------------------------------------------------------------------

    @property
    def cells_count(self):
        """
        Cells count.

        Returns
        -------
        int
            Count of cells.
        """

        return self.a * self.b * self.c

    #-----------------------------------------------------------------------------------------------

    @property
    def weight(self):
        """
        Weight.

        Returns
        -------
        int
            Weight.
        """

        return self.cells_count

    #-----------------------------------------------------------------------------------------------

    def cut(self, p):
        """
        Cut by position.
        We cut only in a-direction.

        Parameters
        ----------
        p : int
            Position.

        Returns
        -------
        (Block, Block)
            Result blocks.
        """

        a, b, c = self.a, self.b, self.c

        assert (p > 0) and (p < a)

        return Block(p, b, c), Block(a - p, b, c)

    #-----------------------------------------------------------------------------------------------

    def cut_half(self):
        """
        Cut half.

        Returns
        -------
        (Block, Block)
            Result blocks.
        """

        return self.cut(self.a // 2)

#===================================================================================================

class Blocks:
    """
    Blocks.
    Blocks are always in descending order.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Constructor.
        """

        self.items = []

    #-----------------------------------------------------------------------------------------------

    def copy(self):
        """
        Get blocks copy.

        Returns
        -------
        Blocks
            Blocks copy.
        """

        bs = Blocks()
        for b in self.items:
            bs.add(b.copy(), is_sort=False)
        bs.sort()

        return bs

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def random(m, alo, ahi, blo, bhi, clo, chi):
        """
        Random blocks set.

        Parameters
        ----------
        m : int
            Blocks count.
        alo : int
            Lo value of side a.
        ahi : int
            Hi value of side a.
        blo : int
            Lo value of side b.
        bhi : int
            Hi value of side b.
        clo : int
            Lo value of side c.
        chi : int
            Hi value of side c.

        Returns
        -------
        Blocks
            Blocks set.
        """

        bs = Blocks()
        for _ in range(m):
            bs.add(Block.random(alo, ahi, blo, bhi, clo, chi), is_sort=False)
        bs.sort()

        return bs

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Bs(count={self.count},cells_count={self.cells_count})'

    #-----------------------------------------------------------------------------------------------

    @property
    def count(self):
        """
        Count of blocks.

        Returns
        -------
        int
            Blocks count.
        """

        return len(self.items)

    #-----------------------------------------------------------------------------------------------

    @property
    def is_empty(self):
        """
        Check if blocks list is empty.

        Returns
        -------
        bool
            True - if blocks list is empty,
            False - otherwise.
        """

        return self.count == 0

    #-----------------------------------------------------------------------------------------------

    def sort(self):
        """
        Sort.
        """

        self.items.sort(key=lambda b: b.cells_count, reverse=True)

    #-----------------------------------------------------------------------------------------------

    def add(self, b, is_sort=True):
        """
        Add block.

        Parameters
        ----------
        b : Block
            Block.
        is_sort : bool
            It is needed to sort after adding.
        """

        self.items.append(b)

        if is_sort:
            self.sort()

    #-----------------------------------------------------------------------------------------------

    def add_blocks(self, bs):
        """
        Add blocks.

        Parameters
        ----------
        bs : [Block]
            List of blocks.
        """

        for b in bs:
            self.add(b, is_sort=False)
        self.sort()

    #-----------------------------------------------------------------------------------------------

    @property
    def cells_count(self):
        """
        Cells count.

        Returns
        -------
        int
            Cells count.
        """

        return sum([b.cells_count for b in self.items])

    #-----------------------------------------------------------------------------------------------

    @property
    def weight(self):
        """
        Weight.

        Returns
        -------
        int
            Weight.
        """

        return self.cells_count

    #-----------------------------------------------------------------------------------------------

    @property
    def first(self):
        """
        First block.

        Returns
        -------
        Block
            First block.
        """

        return self.items[0]

    #-----------------------------------------------------------------------------------------------

    def pop(self, i):
        """
        Pop block.

        Parameters
        ----------
        i : int
            Index of popped block.

        Returns
        -------
        Block
            Popped block.
        """

        return self.items.pop(i)

    #-----------------------------------------------------------------------------------------------

    def pop_first(self):
        """
        Pop first block and return it.

        Returns
        -------
        Block
            First block.
        """

        return self.pop(0)

    #-----------------------------------------------------------------------------------------------

    def remove(self, b):
        """
        Remove blocks.

        Parameters
        ----------
        b : Block
            Block.
        """

        self.items.remove(b)

    #-----------------------------------------------------------------------------------------------

    def clear(self):
        """
        Clear blocks.
        """

        self.items.clear()

#===================================================================================================

class Partition (Blocks):
    """
    Partition.
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

        return f'P({self.weight}w)'

#===================================================================================================

class Partitions:
    """
    Partitions class.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, k):
        """
        Constructor.

        Parameters
        ----------
        k : int
            Partitions count.
        """

        self.items = []
        for _ in range(k):
            self.items.append(Partition())

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String representation.
        """

        return f'Ps(count={self.count},w={self.weight},D*={self.d_star})'

    #-----------------------------------------------------------------------------------------------

    @property
    def count(self):
        """
        Partitions count.

        Returns
        -------
        int
            Partitions count.
        """

        return len(self.items)

    #-----------------------------------------------------------------------------------------------

    @property
    def weight(self):
        """
        Weight.

        Returns
        -------
        int
            Weight.
        """

        return sum([p.weight for p in self.items])

    #-----------------------------------------------------------------------------------------------

    @property
    def lightest(self):
        """
        Get lightest partition.

        Returns
        -------
        Partition
            The lightest partition.
        """

        p = None

        for pi in self.items:
            if p is None:
                p = pi
            elif pi.weight < p.weight:
                p = pi

        return p

    #-----------------------------------------------------------------------------------------------

    @property
    def d_star(self):
        """
        D*.

        Returns
        -------
        float
            D* value.
        """

        ws = [p.weight for p in self.items]
        mw = max(ws)
        aw = sum(ws) / len(self.items)

        if aw > 0.0:
            return (mw - aw) / aw
        else:
            return 0.0

    #-----------------------------------------------------------------------------------------------

    def clear(self):
        """
        Clear all partitions.
        """

        for p in self.items:
            p.clear()

    #-----------------------------------------------------------------------------------------------

    def transfer_all_blocks(self, bs):
        """
        Pass blocks.

        Parameters
        ----------
        bs : Blocks
            Blocks.
        """

        for p in self.items:
            for b in p.items:
                bs.add(b, is_sort=False)
        bs.sort()
        self.clear()

#===================================================================================================

# Distribution methods:
# - greedy - with no cuts;
# - half_max_block - cut max block by half if D* is not reached.

#---------------------------------------------------------------------------------------------------

def distribute_greedy(bs, ps):
    """
    Distribute blocks with greedy strategy.

    Parameters
    ----------
    bs : Blocks
        Blocks.
    ps : Partitions
        Partitions.

    Returns
    -------
    int
        Cuts count.
    """

    while not bs.is_empty:
        ps.lightest.add(bs.pop_first())

    return 0

#---------------------------------------------------------------------------------------------------

def distribute_half_max_block(bs, ps, d_star):
    """
    Distribute block using half max block strategy.

    Parameters
    ----------
    bs : Blocks
        Blocks.
    ps : Partitions
        Partitions.
    d_star : float
        Target D* value.

    Returns
    -------
    int
        Cuts count.
    """

    cc = 0

    while True:
        distribute_greedy(bs, ps)
        if ps.d_star <= d_star:
            break
        ps.transfer_all_blocks(bs)
        bs.add_blocks(bs.pop_first().cut_half())
        cc = cc + 1

    return cc

#---------------------------------------------------------------------------------------------------

def distribute_min_blocks_cuts_find_next(bs, ps, t, mrg, mb,
                                         is_full_under, is_part_under,
                                         is_full_above, is_part_above):
    """
    Find next possible block (or block part) placement into some partition.

    Parameters
    ----------
    bs : Blocks
        Blocks.
    ps : Partitions
        Partitions.
    t : float
        Target weight.
    mrg : int
        Margin of cut.
    mb : int
        Minimal block size.
    is_full_under : bool
        Flag for check place full block under.
    is_part_under
        Flag for check place part block under.
    is_full_above : bool
        Flag for check place full block above.
    is_part_above : bool
        Flag for check place part block above.

    Returns
    -------
    (char, Block, Partition, int):
        Kind of placement, block, partition, cut position.
        Kind of placement may be one of the following values:
            'U' - full under,
            'u' - part under,
            'A' - full above,
            'a' - part above,
            'n' - no placement.
    ('n', None, None, 0):
        If placement can not be found.
    """

    # Initial values.
    k, b, p, c, d = 'n', None, None, 0, math.inf

    # Full under.
    if is_full_under:
        for pi in ps.items:
            pw = pi.weight
            for bi in bs.items:
                bw = bi.weight
                if pw + bw <= t:
                    di = t - (pw + bw)
                    if di < d:
                        k, b, p, c, d = 'U', bi, pi, 0, di

    # Part under.
    if (k == 'n') and is_part_under:
        for pi in ps.items:
            pw = pi.weight
            for bi in bs.items:
                bw = bi.weight
                for ci in range(mrg, bi.a - mrg + 1):
                    bpw = ci * bi.b * bi.c
                    bpw1 = bw - bpw
                    if (bpw >= mb) and (bpw1 >= mb):
                        if pw + bpw <= t:
                            di = t - (pw + bpw)
                            if di < d:
                                k, b, p, c, d = 'u', bi, pi, ci, di

    # Part above.
    if (k == 'n') and is_part_above:
        for pi in ps.items:
            pw = pi.weight
            for bi in bs.items:
                bw = bi.weight
                for ci in range(mrg, bi.a - mrg + 1):
                    bpw = ci * bi.b * bi.c
                    bpw1 = bw - bpw
                    if (bpw >= mb) and (bpw1 >= mb):
                        if pw + bpw > t:
                            di = (pw + bpw) - t
                            if di < d:
                                k, b, p, c, d = 'a', bi, pi, ci, di

    # Full above.
    if (k == 'n') and is_full_above:
        for pi in ps.items:
            pw = pi.weight
            for bi in bs.items:
                bw = bi.weight
                if pw + bw > t:
                    di = (pw + bw) - t
                    if di < d:
                        k, b, p, c, d = 'A', bi, pi, 0, di

    return k, b, p, c

def distribute_min_blocks_cuts(bs, ps, mrg, mb,
                               is_full_under, is_part_under,
                               is_full_above, is_part_above):
    """
    Distribution with blocks cuts minimization.

    Parameters
    ----------
    bs : Blocks
        Blocks.
    ps : Partitions.
        Partitions.
    mrg : int
        Limit margin for result blocks.
    mb : int
        Limit minimal block size.
    is_full_under : bool
        Flag for check place full block under.
    is_part_under
        Flag for check place part block under.
    is_full_above : bool
        Flag for check place full block above.
    is_part_above : bool
        Flag for check place part block above.

    Variants of flags:
    +---------------+---------------+---------------+---------------+---------
    | is_full_under | is_part_under | is_full_above | is_part_above | comment
    +---------------+---------------+---------------+---------------+---------
    | False         | False         | False         | False         | NA - no variants at all
    | False         | False         | False         | True          | NA - bad for small blocks
    | False         | False         | True          | False         | NA - bad for small blocks
    | False         | False         | True          | True          | NA - bad for small blocks
    | False         | True          | False         | False         | NA - bad for big blocks
    | False         | True          | False         | True          | NA - bad if cuts are not available due margin or min block size
    | False         | True          | True          | False         | NA - bad if cuts are not available due margin or min block size
    | False         | True          | True          | True          | NA - bad if cuts are not available due margin or min block size
    | True          | False         | False         | False         | NA - bad for big blocks
    | True          | False         | False         | True          | NA - bad if cuts are not available due margin or min block size
    | True          | False         | True          | False         | BAD - because there is no cuts
    | True          | False         | True          | True          | YES
    | True          | True          | False         | False         | NA - bad for big blocks
    | True          | True          | False         | True          | NA - bad if cuts are not available due margin or min block size
    | True          | True          | True          | False         | YES
    | True          | True          | True          | True          | YES
    +---------------+---------------+---------------+---------------+---------

    Returns
    -------
    int
        Blocks cuts count.
    """

    cc = 0
    t = bs.cells_count / ps.count

    while not bs.is_empty:
        k, b, p, c = distribute_min_blocks_cuts_find_next(bs, ps, t, mrg, mb,
                                                          is_full_under, is_part_under,
                                                          is_full_above, is_part_above)
        if k == 'U':
            p.add(b)
            bs.remove(b)
        elif k == 'u':
            b1, b2 = b.cut(c)
            cc = cc + 1
            p.add(b1, is_sort=False)
            bs.remove(b)
            bs.add(b2)
        elif k == 'A':
            p.add(b)
            bs.remove(b)
            t = p.weight
        elif k == 'a':
            b1, b2 = b.cut(c)
            cc = cc + 1
            p.add(b1, is_sort=False)
            bs.remove(b)
            bs.add(b2)
            t = p.weight
        else:
            print('bs:', bs, bs.items)
            print('ps:', ps, ps.items)
            print('t:', t, '\nmrg:', mrg, '\nmb:', mb)
            print('is_full_under:', is_full_under, '\nis_part_under:', is_part_under)
            print('is_full_above:', is_full_above, '\nis_part_above:', is_part_above)
            assert False

    return cc

#===================================================================================================

def test():
    """
    Tests.
    """

    if MIN_CUTS_FOR_EXTRACT_PART_ON:
        # Minimal cuts count to extract from block n x 1 x 1 part of size t.
        # 1d
        assert math.isinf(min_cuts_for_extract_part_1d(5, 6))
        assert min_cuts_for_extract_part_1d(5, 0) == 0
        assert min_cuts_for_extract_part_1d(5, 5) == 0
        assert min_cuts_for_extract_part_1d(5, 3) == 1
        # 2d
        assert math.isinf(min_cuts_for_extract_part_2d(3, 3, 10))
        assert min_cuts_for_extract_part_2d(3, 3, 0) == 0
        assert min_cuts_for_extract_part_2d(3, 3, 9) == 0
        assert min_cuts_for_extract_part_2d(3, 3, 6) == 1
        assert min_cuts_for_extract_part_2d(3, 3, 4) == 2
        assert min_cuts_for_extract_part_2d(3, 3, 5) == 2
        # 3d
        assert math.isinf(min_cuts_for_extract_part_3d(3, 3, 3, 28))
        assert min_cuts_for_extract_part_3d(3, 3, 3, 0) == 0
        assert min_cuts_for_extract_part_3d(3, 3, 3, 27) == 0
        assert min_cuts_for_extract_part_3d(3, 3, 3, 9) == 1
        assert min_cuts_for_extract_part_3d(5, 5, 5, 50) == 1
        assert min_cuts_for_extract_part_3d(6, 6, 6, 53) == 5
        assert min_cuts_for_extract_part_3d(6, 6, 6, 59) == 5
        assert min_cuts_for_extract_part_3d(6, 6, 6, 89) == 5

#---------------------------------------------------------------------------------------------------

def test_distribute_blocks_step(m, k, lo, hi):
    """
    Blocks distribution test.

    Parameters
    ----------
    m : int
        Blocks count.
    k : int
        Processes count.
    lo : int
        Lo bound of block size.
    hi : int
        Hi bound of block size.

    Returns
    -------
    [gstar, 1011star, 1110star, 1111star, hstar], [gcc, 1011cc, 1110cc, 1111cc, hcc]
        gstar - D* for greedy distribution,
        1011star - D* for 1011 algorithm,
        1110star - D* for 1110 algorithm,
        1111star - D* for 1111 algorithm,
        hstar - D* for half max block allgorithm,
        gcc - cuts count for greedy distribution,
        1011cc - cuts count for 1011 algorithm,
        1110cc - cuts count for 1110 algorithm,
        1111cc - cuts count for 1111 algorithm,
        hcc - cuts count for half max block allgorithm.
    """

    random.seed()

    # Create blocks.
    print('Blocks')
    bs_g = Blocks.random(m, lo, hi, lo, hi, lo, hi)
    bs_m_1011 = bs_g.copy()
    bs_m_1110 = bs_g.copy()
    bs_m_1111 = bs_g.copy()
    bs_h = bs_g.copy()
    print(bs_g.items)
    print('bs_g', bs_g)

    # Distribution parameters.
    margin = 3
    min_block = 64
    ps_g = Partitions(k)
    ps_m_1011 = Partitions(k)
    ps_m_1110 = Partitions(k)
    ps_m_1111 = Partitions(k)
    ps_h = Partitions(k)

    # Distribute blocks.
    cc_g = distribute_greedy(bs_g, ps_g)
    cc_m_1011 = distribute_min_blocks_cuts(bs_m_1011, ps_m_1011, margin, min_block,
                                           True,  False,
                                           True,  True)
    cc_m_1110 = distribute_min_blocks_cuts(bs_m_1110, ps_m_1110, margin, min_block,
                                           True,  True,
                                           True,  False)
    cc_m_1111 = distribute_min_blocks_cuts(bs_m_1111, ps_m_1111, margin, min_block,
                                           True,  True,
                                           True,  True)
    d_star = min([ps_m_1011.d_star, ps_m_1110.d_star, ps_m_1111.d_star])
    cc_h = distribute_half_max_block(bs_h, ps_h, d_star)

    # Print result.
    print('Partitions:')
    print('Greedy     :', ps_g, cc_g)
    print('MinCut 1011:', ps_m_1011, cc_m_1011)
    print('MinCut 1110:', ps_m_1110, cc_m_1110)
    print('MinCut 1111:', ps_m_1111, cc_m_1111)
    print('HalfBl     :', ps_h, cc_h)

    return [ps_g.d_star, ps_m_1011.d_star, ps_m_1110.d_star, ps_m_1111.d_star, ps_h.d_star], \
           [cc_g, cc_m_1011, cc_m_1110, cc_m_1111, cc_h]

#---------------------------------------------------------------------------------------------------

def test_distribute_blocks_steps(runs, m, k, lo, hi):
    """
    Blocks distribution test.

    Parameters
    ----------
    runs : int
        Runs count.
    m : int
        Blocks count.
    k : int
        Processes count.
    lo : int
        Lo bound of block size.
    hi : int
        Hi bound of block size.
    """

    s_g, s_m_1011, s_m_1110, s_m_1111, s_h = [], [], [], [], []
    c_g, c_m_1011, c_m_1110, c_m_1111, c_h = [], [], [], [], []
    for i in range(runs):
        print('STEP :', i)
        r = test_distribute_blocks_step(m, k, lo, hi)
        s_g.append(r[0][0])
        s_m_1011.append(r[0][1])
        s_m_1110.append(r[0][2])
        s_m_1111.append(r[0][3])
        s_h.append(r[0][4])
        c_g.append(r[1][0])
        c_m_1011.append(r[1][1])
        c_m_1110.append(r[1][2])
        c_m_1111.append(r[1][3])
        c_h.append(r[1][4])
    print('s_g =', s_g)
    print('s_m_1011 =', s_m_1011)
    print('s_m_1110 =', s_m_1110)
    print('s_m_1111 =', s_m_1111)
    print('s_h =', s_h)
    print('c_g =', c_g)
    print('c_m_1011 =', c_m_1011)
    print('c_m_1110 =', c_m_1110)
    print('c_m_1111 =', c_m_1111)
    print('c_h =', c_h)

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    test()

    if MIN_CUTS_FOR_EXTRACT_PART_ON:
        # Draw plot for P(n, m, k, t).
        plt.figure(figsize=(10, 6), dpi=100)
        n = 6
        xs = range(n**3 + 1)
        plt.plot(xs, [min_cuts_for_extract_part_3d(n, n, n, x) for x in xs],
                 marker='o', linewidth=2.0)
        plt.xlabel(r'Размер выделяемой части $t$')
        plt.ylabel(r'Минимальное количество разрезов $P_{6 \times 6 \times 6}^{t}$')
        plt.grid(True)
        plt.show()

    test_distribute_blocks_steps(10, 10, 10, 20, 30)

#===================================================================================================
