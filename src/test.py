"""
Test.
"""

import mesh_su_self_intersections_elimination

#===================================================================================================

# Test self-intersections elimination.
mesh_name = '../data/meshes/tu/tu1_078_int_closed'
mesh_su_self_intersections_elimination.sie(f'{mesh_name}.dat',
                                           f'{mesh_name}_out.dat',
                                           1000000)

#===================================================================================================
