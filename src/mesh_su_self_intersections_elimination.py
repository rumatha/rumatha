#!/usr/bin/env python3
"""
Self-intersections elimination from mesh.
"""

import argparse
import mesh_su
import geom3d_rat
import time

#===================================================================================================

def sie(in_mesh, out_mesh, denom):
    """
    Delete self-intersections.

    Parameters
    ----------
    in_mesh : str
        In mesh file.
    out_mesh : str
        Out mesh file.
    denom : int
        Denominator for rational numbers.
    """

    start = time.time()
    mesh = mesh_su.Mesh(in_mesh)

    # Delete self-intersections.
    m = mesh.delete_self_intersections_rat(denom=denom, is_log=True)
    m.store(out_mesh)

    # Print stat.
    geom3d_rat.print_statistics()
    print(f'total time : {(time.time() - start):.4}')

#===================================================================================================

if __name__ == '__main__':

    # Arguments parser.
    parser = argparse.ArgumentParser(prog='mesh_su_delete_self_intersections',
                                     description='Delete self-intersections from mesh.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_mesh', help='Input mesh file.')
    parser.add_argument('out_mesh', help='Output mesh file.')
    parser.add_argument('denom', help='Denominator for rational numbers.', type=int)
    args = parser.parse_args()

    # Print parameters.
    print('Input data:')
    print(f'\tin_mesh  : {args.in_mesh}')
    print(f'\tout_mesh : {args.out_mesh}')
    print(f'\tdenom    : {args.denom}')

    # Run.
    sie(args.in_mesh, args.out_mesh, args.denom)

#===================================================================================================
