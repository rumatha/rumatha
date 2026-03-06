"""
Mesh - Surface Unstructured.
"""

import numpy as np
import geom3d_rat
from fractions import Fraction as Fr
import time

#===================================================================================================

# Valuable digits in coordinates.
NODE_COORDINATES_VALUABLE_DIGITS_COUNT = 10

# String of export.
EXPORT_FORMAT_STRING = '{0:.18e}'

#===================================================================================================

class Node:
    """
    Node - container for coordinates.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, p):
        """
        Initialization.
        Node may appear only as point holder.

        Parameters
        ----------
        p : np.array
            Point coordinates.
        """

        # Global identifier.
        self.glo_id = -1

        self.p = p

        # Incident objects.
        self.edges = []
        self.faces = []

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String.
        """

        return f'Node {self.glo_id} ({self.p})'

    #-----------------------------------------------------------------------------------------------

    def rounded_coordinates(self):
        """
        Tuple with rounded coordinates.

        Returns
        -------
        tuple
            Rounded coordinates.
        """

        p = self.p
        if isinstance(p, geom3d_rat.Point):
            p = p.get_real_array()

        return tuple(map(lambda x: round(x, NODE_COORDINATES_VALUABLE_DIGITS_COUNT), p))

    #-----------------------------------------------------------------------------------------------

    def is_isolated(self):
        """
        Check if node is isolated.

        Returns
        -------
        True - if node is isolated,
        False - otherwise.
        """

        return len(self.edges) == 0

    #-----------------------------------------------------------------------------------------------

    def is_border(self):
        """
        Check if node is border.

        Returns
        -------
        bool
            True - if node is border,
            False - otherwise.
        """

        return any(e.is_border() for e in self.edges)

    #-----------------------------------------------------------------------------------------------

    def neighbour(self, e):
        """
        Get neighbour by edge.

        Parameters
        ----------
        e : Edge
            Edge.

        Returns
        -------
        Node
            Neighbour node or None.
        """

        if self == e.nodes[0]:
            return e.nodes[1]
        elif self == e.nodes[1]:
            return e.nodes[0]
        else:
            return None

    #-----------------------------------------------------------------------------------------------

    def neighbourhood(self):
        """
        Get heighbourhood.

        Returns
        -------
        [Node]
            List of neighbour nodes.
        """

        return [self.neighbour(e) for e in self.edges]

#===================================================================================================

class Edge:
    """
    Edge - border between two faces
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Initialization.
        """

        # Global identifier.
        self.glo_id = -1

        # Incident objects.
        self.faces = []
        self.nodes = []

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String.
        """

        id0, id1 = self.nodes[0].glo_id, self.nodes[1].glo_id

        return f'Edge {self.glo_id} ({id0} - {id1})'

    #-----------------------------------------------------------------------------------------------

    def is_faces_free(self):
        """
        Check if edge is without incident faces.

        Returns
        -------
        True - if edge is without faces,
        False - otherwise.
        """

        return len(self.faces) == 0

    #-----------------------------------------------------------------------------------------------

    def is_border(self):
        """
        Check if edge is border.

        Returns
        -------
        bool
            True - if edge is border,
            False - otherwise.
        """

        return len(self.faces) == 1

    #-----------------------------------------------------------------------------------------------

    def points(self):
        """
        Get points.

        Returns
        -------
        (Point, Point)
            Points.
        """
        return self.nodes[0].p, self.nodes[1].p

    #-----------------------------------------------------------------------------------------------

    def flip_nodes(self):
        """
        Flip nodes.
        """

        self.nodes[0], self.nodes[1] = self.nodes[1], self.nodes[0]

#===================================================================================================

class Face:
    """
    Face - container for physical data.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self):
        """
        Initialization.
        """

        # Global identifier.
        self.glo_id = -1

        # Data.
        self.data = dict()

        # Incident objects.
        self.nodes = []
        self.edges = []

        # Link to zone.
        self.zone = None

    #-----------------------------------------------------------------------------------------------

    def set_data(self, variables, values):
        """
        Set data.

        Parameters
        ----------
        variables : [str]
            Variables names.
        values : [object]
            List of values.

        """

        self.data = dict(zip(variables, values))

    #-----------------------------------------------------------------------------------------------

    def copy_data_from(self, f):
        """
        Copy data from another face.

        Parameters
        ----------
        f : Face
            Another face.
        """

        self.set_data(f.data.keys(), f.data.values())

    #-----------------------------------------------------------------------------------------------

    def __repr__(self):
        """
        String representation.

        Returns
        -------
        str
            String.
        """

        id0, id1, id2 = self.nodes[0].glo_id, self.nodes[1].glo_id, self.nodes[2].glo_id

        return f'Face {self.glo_id} ({id0}, {id1}, {id2})'

    #-----------------------------------------------------------------------------------------------

    def __getitem__(self, item):
        """
        Get face data element.

        Parameters
        ----------
        item : str
            Name of data element.

        Returns
        -------
        value
            Value of data element.
        """

        return self.data.get(item, 0.0)

    #-----------------------------------------------------------------------------------------------

    def __setitem__(self, key, value):
        """
        Set data element.

        Parameters
        ----------
        key : str
            Name of data element.
        value
            Value of data element.
        """

        self.data[key] = value

    #-----------------------------------------------------------------------------------------------

    def points(self):
        """
        Get points.

        Returns
        -------
        tuple
            Points.
        """

        return self.nodes[0].p, self.nodes[1].p, self.nodes[2].p

    #-----------------------------------------------------------------------------------------------

    def outer_normal(self):
        """
        Outer normal.

        Returns
        -------
        np.array
            Outer normal.
        """

        return np.cross(self.nodes[1].p - self.nodes[0].p, self.nodes[2].p - self.nodes[0].p)

    #-----------------------------------------------------------------------------------------------

    def neighbour(self, e):
        """
        Get neighbour by edge.

        Parameters
        ----------
        e : Edge
            Edge.

        Returns
        -------
        Face | None
            Neighbour face or None.
        """

        assert len(e.faces) == 2

        if self == e.faces[0]:
            return e.faces[1]
        elif self == e.faces[1]:
            return e.faces[0]
        else:
            return None

    #-----------------------------------------------------------------------------------------------

    def neighbourhood(self):
        """
        Get neighbourhood.

        Returns
        -------
        [Face]
            List of neighbour faces (by all edges).
        """

        nh = []

        for e in self.edges:
            for f in e.faces:
                if f != self:
                    nh.append(f)

        return nh

    #-----------------------------------------------------------------------------------------------

    def third_node(self, e):
        """
        Get third node, not incident to edge.

        Parameters
        ----------
        e : Edge
            Edge.

        Returns
        -------
        Node
            Third node.
        """

        a, b = e.nodes[0], e.nodes[1]

        for n in self.nodes:
            if (n != a) and (n != b):
                return n

        assert False

#===================================================================================================

class Zone:
    """
    Zone - set of faces.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, name):
        """
        Initialization.

        Parameters
        ----------
        name : str
            Name of zone.
        """

        self.name = name
        self.nodes = []
        self.faces = []

    #-----------------------------------------------------------------------------------------------

    @staticmethod
    def objects_slice_str(fun, obs):
        """
        String, that contains data slice for some objects.
        Formatting is defined here.

        Parameters
        ----------
        fun
            Function for data extracting.
        obs
            Objects list.

        Returns
        -------
        str
            String with data slice.
        """

        return ' '.join(map(lambda ob: EXPORT_FORMAT_STRING.format(fun(ob)), obs))

    #-----------------------------------------------------------------------------------------------

    def nodes_coordinate_slice_str(self, i):
        """
        String, that contains i-th coordinate of all nodes.

        Parameters
        ----------
        i : int
            Coordinate index.

        Returns
        -------
        str
            String with coordinate slice.

        """

        return Zone.objects_slice_str(lambda n: n.p[i], self.nodes)

    #-----------------------------------------------------------------------------------------------

    def faces_data_element_slice_str(self, e):
        """
        String, that contains data element for all faces.

        Parameters
        ----------
        e : str
            Name of data element.

        Returns
        -------
        str
            String with data element slice.
        """

        return Zone.objects_slice_str(lambda f: f[e], self.faces)

#===================================================================================================

class Mesh:
    """
    Mesh - consists of surface triangle faces.
    """

    #-----------------------------------------------------------------------------------------------

    def __init__(self, filename=None, is_merge_nodes=True):
        """
        Initialization.

        Parameters
        ----------
        filename : str
            File for load.
        is_merge_nodes : bool
            Is merge nodes flag.
        """

        # Comment and title - save for store.
        self.comment = ''
        self.title = ''

        # Set empty sets of nodes, faces, zones.
        self.zones = []
        self.nodes = []
        self.faces = []
        self.edges = []

        # Rounded coordinates bag.
        self.rounded_coordinates_bag = set()

        # Load.
        if not filename is None:
            self.load(filename, is_merge_nodes=is_merge_nodes)

    #-----------------------------------------------------------------------------------------------

    def clear(self):
        """
        Clear all.
        """

        self.comment = ''
        self.title = ''
        self.nodes.clear()
        self.edges.clear()
        self.faces.clear()
        self.zones.clear()
        self.rounded_coordinates_bag.clear()

    #-----------------------------------------------------------------------------------------------

    def print(self,
              print_edges_with_incident_faces=False,
              print_faces_neighbourhood=False):
        """
        Print information.

        Parameters
        ---------
        print_edges_with_incident_faces : bool
            Flag for print edges with incident faces.
        print_faces_neighbourhood : bool
            Flag for print faces with neighbourhood.
        """

        print('[MESH]')
        print(f'Nodes ({len(self.nodes)}):\n  ', self.nodes)
        print(f'Edges ({len(self.edges)}):\n  ', self.edges)
        print(f'Faces ({len(self.faces)}):\n  ', self.faces)

        if print_edges_with_incident_faces:
            print('[EDGES WITH INCIDENT FACES]')
            for e in self.edges:
                print(f'{e} --- {len(e.faces)} / {e.faces}')

        if print_faces_neighbourhood:
            print('[FACES WITH NEIGHBOURHOOD]')
            for f in self.faces:
                nh = f.neighbourhood()
                print(f'{f} --- {len(nh)} / {nh}')

    #-----------------------------------------------------------------------------------------------

    def find_near_node(self, node):
        """
        Try to find node near to given node.

        Parameters
        ----------
        node : Node
            Given node.

        Returns
        -------
        Node or None
            If node is found, return it, otherwise return None.
        """

        rc = node.rounded_coordinates()

        # Try to find in bag.
        if rc not in self.rounded_coordinates_bag:
            return None

        # Node rounded coordinates is in bag, find it.
        for n in self.nodes:
            if rc == n.rounded_coordinates():
                return n

        raise Exception('Internal error')

    #-----------------------------------------------------------------------------------------------

    def find_edge(self, a, b, except_edge=None):
        """
        Find edge with two nodes.

        Parameters
        ----------
        a : Node
            First node.
        b : Node.
            Second node.
        except_edge: Edge
            parameter for searching double edges

        Returns
        -------
        Edge | None
            Found edge or None.
        """

        for e in a.edges:
            if a.neighbour(e) == b and e != except_edge:
                return e

        # Not found.
        return None

    #-----------------------------------------------------------------------------------------------

    def find_face(self, a, b, c):
        """
        Find face with given nodes.

        Parameters
        ----------
        a : Node
            First node.
        b : Node
            Second node.
        c : Node
            Third node.

        Returns
        -------
        Face | None
        """

        ids = sorted([a.glo_id, b.glo_id, c.glo_id])

        for f in a.faces:
            lids = sorted([n.glo_id for n in f.nodes])
            if ids == lids:
                return f

        return None

    #-----------------------------------------------------------------------------------------------

    def max_node_glo_id(self):
        """
        Get maximum node global id
        (id of the last node).

        Returns
        -------
        int
            Maximum node global id,
            or -1, if there is no nodes.
        """

        if self.nodes:
            return self.nodes[-1].glo_id
        else:
            return -1

    #-----------------------------------------------------------------------------------------------

    def max_edge_glo_id(self):
        """
        Get maximum edge global id
        (id of the last edge).

        Returns
        -------
        int
            Maximum edge global id,
            or -1, if there is no edges.
        """

        if self.edges:
            return self.edges[-1].glo_id
        else:
            return -1

    #-----------------------------------------------------------------------------------------------

    def max_face_glo_id(self):
        """
        Get maximum face global id
        (id of the last face).

        Returns
        -------
        int
            Maximum face global if,
            or -1, if there is no faces.
        """

        if self.faces:
            return self.faces[-1].glo_id
        else:
            return -1

    #-----------------------------------------------------------------------------------------------

    def add_zone(self, name):
        """
        Add zone.

        Parameters
        ----------
        name : str
            Name of zone.

        Returns
        -------
        Zone
            Added zone.
        """

        z = Zone(name)
        self.zones.append(z)

        return z

    #-----------------------------------------------------------------------------------------------

    def add_node(self, p, zone, is_merge_nodes=True):
        """
        Add node to mesh.
        This is only way to add node into mesh.

        Parameters
        ----------
        p : Point
            Point.
        zone : Zone
            Zone to add node to.
        is_merge_nodes : bool
            Is merge nodes flag.

        Returns
        -------
        Node
            Added node
            (it may be new node or found near node).
        """

        n = Node(p)

        if is_merge_nodes:
            found_node = self.find_near_node(n)
        else:
            found_node = None

        if found_node is None:
            max_glo_id = self.max_node_glo_id()
            n.glo_id = max_glo_id + 1
            self.nodes.append(n)
            self.rounded_coordinates_bag.add(n.rounded_coordinates())
            node_to_zone = n
        else:
            node_to_zone = found_node

        zone.nodes.append(node_to_zone)

        return node_to_zone

    #-----------------------------------------------------------------------------------------------

    def add_edge(self, a, b):
        """
        Add edge or return already existing one.

        Parameters
        ----------
        a : Node
            First node.
        b : Node
            Second node.

        Returns
        -------
        Edge
            Edge (found or new).
        """

        e = self.find_edge(a, b)

        if e is None:
            e = Edge()
            max_glo_id = self.max_edge_glo_id()
            e.glo_id = max_glo_id + 1
            self.edges.append(e)
            self.links([(a, e), (b, e)])

        return e

    #-----------------------------------------------------------------------------------------------

    def add_face(self, a, b, c, zone):
        """
        Add face to mesh.

        Parameters
        ----------
        a : Node
            First node.
        b : Node
            Second node.
        c : Node
            Third node.
        zone : Zone
            Zone to add to.
        """

        f = self.find_face(a, b, c)

        if not f is None:
            return f

        f = Face()
        max_glo_id = self.max_face_glo_id()
        f.glo_id = max_glo_id + 1
        self.faces.append(f)
        zone.faces.append(f)
        f.zone = zone
        ab, bc, ac = self.add_edge(a, b), self.add_edge(b, c), self.add_edge(a, c)
        self.links([(a, f), (b, f), (c, f), (ab, f), (bc, f), (ac, f)])

        return f

    #-----------------------------------------------------------------------------------------------

    def add_triangle(self, z, t):
        """
        Add triangle to zone.

        Parameters
        ----------
        z : Zone
            Zone.
        t : Triangle
            Triangle.
        """

        a = self.add_node(t.A, z, True)
        b = self.add_node(t.B, z, True)
        c = self.add_node(t.C, z, True)
        self.add_face(a, b, c, z)

    #-----------------------------------------------------------------------------------------------

    def add_triangles(self, z, ts):
        """
        Add triangles to zone.

        Parameters
        ----------
        z : Zone
            Zone.
        ts : Triangles
            Triangles.
        """

        for t in ts:
            self.add_triangle(z, t)

    #-----------------------------------------------------------------------------------------------

    def link(self, obj1, obj2):
        """
        Link two objects.
        Objects that can be linked:
          - node - edge
          - node - face
          - edge - face

        Parameters
        ----------
        obj1 : Node | Edge
            First object.
        obj2 : Edge | Face
            Second object.
        """

        if isinstance(obj1, Node):
            if isinstance(obj2, Edge):
                obj1.edges.append(obj2)
                obj2.nodes.append(obj1)
            elif isinstance(obj2, Face):
                obj1.faces.append(obj2)
                obj2.nodes.append(obj1)
            else:
                raise Exception(f'msu.Mesh : wrong object type in link ({obj2})')
        elif isinstance(obj1, Edge):
            if isinstance(obj2, Face):
                obj1.faces.append(obj2)
                obj2.edges.append(obj1)
            else:
                raise Exception(f'msu.Mesh : wrong object type in link ({obj2})')
        else:
            raise Exception(f'msu.Mesh : wrong object type in link ({obj1})')

    #-----------------------------------------------------------------------------------------------

    def links(self, li):
        """
        Multiple links.

        Parameters
        ----------
        li : [(object, object]
            List of pairs for link.
        """

        for obj1, obj2 in li:
            self.link(obj1, obj2)

    #-----------------------------------------------------------------------------------------------

    def unlink(self, obj1, obj2):
        """
        Unlink two objects.

        Parameters
        ----------
        obj1 : Node | Edge
            First object.
        obj2 : Edge | Face
            Second object.
        """

        if isinstance(obj1, Node):
            if isinstance(obj2, Edge):
                obj1.edges.remove(obj2)
                obj2.nodes.remove(obj1)
            elif isinstance(obj2, Face):
                obj1.faces.remove(obj2)
                obj2.nodes.remove(obj1)
            else:
                raise Exception(f'msu.Mesh : wrong object type in unlink ({obj2})')
        elif isinstance(obj1, Edge):
            if isinstance(obj2, Face):
                obj1.faces.remove(obj2)
                obj2.edges.remove(obj1)
            else:
                raise Exception(f'msu.Mesh : wrong object type in unlink ({obj2})')
        else:
            raise Exception(f'msu.Mesh : wrong object type in unlink ({obj1})')

    #-----------------------------------------------------------------------------------------------

    def unlinks(self, li):
        """
        Multiple unlink.

        Parameters
        ----------
        li : [(object, object)]
            List  of pairs for unlink.
        """

        for obj1, obj2 in li:
            self.unlink(obj1, obj2)

    #-----------------------------------------------------------------------------------------------

    def create_edges(self):
        """
        Delete all edges and create them.
        """

        # Delete all edges manually.
        # After this action the mesh is not consistent.
        for n in self.nodes:
            n.edges = []
        for f in self.faces:
            f.edges = []
        self.edges = []

        # Construct edges.
        for f in self.faces:
            a, b, c = f.nodes[0], f.nodes[1], f.nodes[2]
            for first, second in [(a, b), (b, c), (a, c)]:
                e = self.add_edge(first, second)
                self.link(e, f)

    #-----------------------------------------------------------------------------------------------

    def load(self, filename, is_merge_nodes=True):
        """
        Load mesh.

        Parameters
        ----------
        filename : str
            Name of file.
        is_merge_nodes : bool
            Is merge nodes flag.
        """

        variables = []
        face_variables = []
        face_variables_count = 0

        # Clear all objects of the grid.
        self.clear()

        # Open file and try to load it line by line.
        with open(filename, 'r') as f:
            line = f.readline()
            while line:

                if line[0] == '#':

                    # Comment, save it.
                    self.comment = line[1:-1]

                elif 'TITLE=' in line:

                    # Title, save it.
                    self.title = line.split('=')[1][1:-2]

                elif 'VARIABLES=' in line:

                    # Variables.
                    variables_str = line.split('=')[1][:-1]
                    variables = variables_str.replace('"', '').replace(',', '').split()
                    face_variables = variables[3:]
                    face_variables_count = len(face_variables)

                elif 'ZONE T=' in line:

                    # New zone.
                    zone_name = line.split('=')[1][1:-2]
                    zone = self.add_zone(zone_name)

                    # Read count of nodes and faces to read.
                    nodes_line = f.readline()
                    faces_line = f.readline()
                    packing_line = f.readline()
                    zonetype_line = f.readline()
                    varlocation_line = f.readline()
                    if 'NODES=' not in nodes_line:
                        raise Exception('Wrong nodes line ({0}).'.format(nodes_line))
                    if 'ELEMENTS=' not in faces_line:
                        raise Exception('Wrong faces line ({0}).'.format(faces_line))
                    if 'DATAPACKING=BLOCK' != packing_line[:-1]:
                        raise Exception('Wrong packing line ({0}).'.format(packing_line))
                    if 'ZONETYPE=FETRIANGLE' != zonetype_line[:-1]:
                        raise Exception('Wrong zonetype line ({0}).'.format(zonetype_line))
                    right_varlocation_line = 'VARLOCATION=' \
                                             '([4-{0}]=CELLCENTERED)'.format(len(variables))
                    if right_varlocation_line != varlocation_line[:-1]:
                        raise Exception('Wrong varlocation line ({0}). '
                                        'Right value is {1}'.format(varlocation_line,
                                                                    right_varlocation_line))
                    nodes_to_read = int(nodes_line.split('=')[1][:-1])
                    faces_to_read = int(faces_line.split('=')[1][:-1])

                    # Read data for nodes.
                    c = []
                    for i in range(3):
                        line = f.readline()
                        c.append([float(xi) for xi in line.split()])
                    for i in range(nodes_to_read):
                        self.add_node(np.array([c[0][i], c[1][i], c[2][i]]), zone, is_merge_nodes=is_merge_nodes)

                    # Read data for faces.
                    d = []
                    for i in range(face_variables_count):
                        line = f.readline()
                        d.append([float(xi) for xi in line.split()])
                    values = [[d[j][i] for j in range(face_variables_count)] for i in range(faces_to_read)]

                    # Read connectivity lists.
                    for i in range(faces_to_read):
                        line = f.readline()
                        nodes = [zone.nodes[int(ss) - 1] for ss in line.split()]
                        assert len(nodes) == 3
                        face = self.add_face(nodes[0], nodes[1], nodes[2], zone)
                        face.set_data(face_variables, values[i])
                else:
                    raise Exception('Unexpected line : {0}.'.format(line))

                line = f.readline()
            f.close()

        # Create edges.
        self.create_edges()

        # Set identifiers.
        for i, f in enumerate(self.faces):
            f['Id'] = i
        for n in self.nodes:
            if len(n.faces) == 0:
                self.nodes.remove(n)

    #-----------------------------------------------------------------------------------------------

    def convert_coordinates_real_to_rat(self, denom):
        """
        Convert all nodes coordinates from real to rat.

        Parameters
        ----------
        denom : int
            Denominator.
        """

        for node in self.nodes:
            node.p = geom3d_rat.Point.from_real_array(node.p, denom)

    #-----------------------------------------------------------------------------------------------

    def convert_coordinates_rat_to_real(self):
        """
        Convert all nodes coordinates from rat to real.
        """

        for node in self.nodes:
            coords = [node.p.x, node.p.y, node.p.z]
            node.p = np.array([f.numerator / f.denominator for f in coords])

    #-----------------------------------------------------------------------------------------------

    def set_faces_variables(self, fv):
        """
        Delete all faces variables and set new (with values 0.0).

        Parameters
        ----------
        fv : [str]
            List of faces variables.
        """

        self.set_faces_variables = fv

        for f in self.faces:
            f.data.clear()
            for v in self.set_faces_variables:
                f[v] = 0.0

    #-----------------------------------------------------------------------------------------------

    def store(self, filename):
        """
        Store mesh.

        Parameters
        ----------
        filename : str
            Name of file.
        """

        if not self.faces:
            print('Mesh.store: empty mesh, nothing to store.')
            return

        # Save faces glo_id.
        for f in self.faces:
            f['Id'] = f.glo_id

        variables = ['X', 'Y', 'Z'] + list(self.faces[0].data.keys())

        with open(filename, 'w', newline='\n') as f:

            # Store head.
            f.write(f'#{self.comment}\n')
            f.write(f'TITLE="{self.title}"\n')
            f.write('VARIABLES={0}\n'.format(', '.join(['"{0}"'.format(k) for k in variables])))

            # Store zones.
            for zone in self.zones:

                # Do not store empty zones.
                if len(zone.nodes) == 0:
                    continue

                # Store zone head.
                f.write(f'ZONE T="{zone.name}"\n')
                f.write(f'NODES={len(zone.nodes)}\n')
                f.write(f'ELEMENTS={len(zone.faces)}\n')
                f.write('DATAPACKING=BLOCK\n')
                f.write('ZONETYPE=FETRIANGLE\n')
                f.write(f'VARLOCATION=([4-{len(variables)}]=CELLCENTERED)\n')

                # Write first 3 data items (X, Y, Z coordinates).
                for i in range(3):
                    f.write(zone.nodes_coordinate_slice_str(i) + ' \n')

                # Write rest faces data items.
                for e in variables[3:]:
                    f.write(zone.faces_data_element_slice_str(e) + ' \n')

                # Write connectivity lists.
                for face in zone.faces:
                    f.write(' '.join([str(zone.nodes.index(n) + 1) for n in face.nodes]) + '\n')

            f.close()

    #-----------------------------------------------------------------------------------------------
    # Delete elements.
    #-----------------------------------------------------------------------------------------------

    def delete_face(self, f):
        """
        Delete face.

        Parameters
        ----------
        f : Face
            Face to delete.
        """

        # Unlink from nodes.
        while f.nodes:
            self.unlink(f.nodes[0], f)

        # Unlink from edges.
        while f.edges:
            self.unlink(f.edges[0], f)

        # Remove from zones.
        for z in self.zones:
            if f in z.faces:
                z.faces.remove(f)

        # Remove from mesh.
        self.faces.remove(f)

    #-----------------------------------------------------------------------------------------------

    def delete_faces(self, p):
        """
        Delete faces with predicate.

        Parameters
        ----------
        p : lambda
            Predicate for delete face.
        """

        fs = [f for f in self.faces if p(f)]

        for f in fs:
            self.delete_face(f)

    #-----------------------------------------------------------------------------------------------

    def delete_edge(self, e):
        """
        Delete edge.

        Parameters
        ----------
        e : Edge
            Edge to delete.
        """

        # First we must to delete incident faces.
        while e.faces:
            self.delete_face(e.faces[0])

        # Unlink edge from nodes.
        while e.nodes:
            self.unlink(e.nodes[0], e)

        # Remove from mesh.
        if e in self.edges:
            self.edges.remove(e)

    #-----------------------------------------------------------------------------------------------

    def delete_edges(self, p):
        """
        Delete all edges with predicate.

        Parameters
        ----------
        p : lambda
            Predicate for edge delete.
        """

        es = [e for e in self.edges if p(e)]

        for e in es:
            self.delete_edge(e)

    #-----------------------------------------------------------------------------------------------

    def delete_faces_free_edges(self):
        """
        Delete faces free edges.
        """

        self.delete_edges(lambda e: e.is_faces_free())

    #-----------------------------------------------------------------------------------------------

    def delete_node(self, n, delete_isolated=True):
        """
        Delete node.

        Parameters
        ----------
        n : Node
            Node to be deleted.
        delete_isolated : Bool
            delete_isolated nodes or not
        """

        # First we must delete all adjacent edges
        while n.edges:
            self.delete_edge(n.edges[0])

        # Remove node from zones.
        for z in self.zones:
            if n in z.nodes:
                z.nodes.remove(n)

        # Remove node from mesh if it still there.
        if n in self.nodes:
            self.nodes.remove(n)

    #-----------------------------------------------------------------------------------------------

    def delete_nodes(self, p):
        """
        Delete nodes with predicate.

        Parameters
        ----------
        p : lambda
            Predicate for delete.
        """

        ns = [n for n in self.nodes if p(n)]

        for n in ns:
            self.delete_node(n)

    #-----------------------------------------------------------------------------------------------

    def delete_isolated_nodes(self):
        """
        Delete isolated nodes.
        """

        self.delete_nodes(lambda n: n.is_isolated())

    #-----------------------------------------------------------------------------------------------
    # Reduce.
    #-----------------------------------------------------------------------------------------------

    def reduce_edge(self, e):
        """
        Reduce edge.

        Parameters
        ----------
        e : Edge
            Edge.
        """

        # Get all objects needed to process.
        [a, b] = e.nodes
        assert a != b

        # Delete edge e
        self.delete_edge(e)

        # Correct coordinate.
        a.p = 0.5 * (a.p + b.p)

        # For all faces incident to node b create twin for node a.
        for f in b.faces:
            ns = [f.nodes[0], f.nodes[1], f.nodes[2]] # create new list of nodes
            ns[ns.index(b)] = a
            # Force add to keep closed surface.
            r = self.add_face(ns[0], ns[1], ns[2], f.zone)

        # Delete extra node b.
        self.delete_node(b)

        # TODO.
        # Delete bad objects (this is extra code, can be simplified).
        self.delete_faces_free_edges()
        self.delete_isolated_nodes()

    #-----------------------------------------------------------------------------------------------
    # Split.
    #-----------------------------------------------------------------------------------------------

    def parallel_move(self, v):
        """
        Parallel move all nodes.

        Parameters
        ----------
        v : Vector
            Move vector.
        """

        for n in self.nodes:
            n.p += v

    # ----------------------------------------------------------------------------------------------

    def unite_with(self, m):
        """
        Unite with another mesh.

        Parameters
        ----------
        m : Mesh
            Mesh.
        """

        # Correct m zones names.
        for z in m.zones:
            z.name = z.name + ' (unite)'

        # Dumb direct merge (may be incorrect).
        self.zones = self.zones + m.zones
        self.nodes = self.nodes + m.nodes
        self.edges = self.edges + m.edges
        self.faces = self.faces + m.faces

    #-----------------------------------------------------------------------------------------------

    def check_mesh_is_closed(self):
        """
        Check mesh is closed.

        Returns
        -------
        bool
            True - if mesh is closed,
            False - otherwise.
        """

        assert all([len(e.faces) == 2 for e in self.edges])

    #-----------------------------------------------------------------------------------------------

    def delete_self_intersections_rat(self, denom, is_log=False):
        """
        Delete self-intersections using rational coordinates.

        Parameters
        ----------
        denom : int
            Denominator for rational coordinates.
        is_log : bool
            Logger flag.

        Returns
        -------
        Mesh
            Result mesh.
        """

        #
        # First phase. Prepare mesh.
        #

        if is_log:
            print('DSI.Phase.1 : prepare : begin')

        # Convert coordinates to rational.
        self.convert_coordinates_real_to_rat(denom)

        # Create collection of triangles and set of intersection points and segments.
        tis = [(geom3d_rat.Triangle(f.nodes[0].p, f.nodes[1].p, f.nodes[2].p),
                geom3d_rat.PointsAndSegments()) for f in self.faces]
        n = len(tis)

        # Create result mesh with single zone.
        m = Mesh()
        z = m.add_zone('SINGLE ZONE')

        if is_log:
            print('DSI.Phase.1 : prepare : end')

        #
        # Second phase. Find triangles pairs and find all intersection points and segments.
        #

        if is_log:
            print('DSI.Phase.2 : intsec : begin')

        # Process every triangles pair.
        for i in range(n):
            (t1, intsec1) = tis[i]
            for j in range(i + 1, n):
                (t2, intsec2) = tis[j]

                # Find intersection.
                r = geom3d_rat.Intersection.triangle_triangle(t1, t2)

                # Analyze types of intersection.
                if r is None:
                    pass
                elif isinstance(r, geom3d_rat.Point):
                    if not r.is_triangle_vertex(t1):
                        intsec1.add_unique_point(r)
                    if not r.is_triangle_vertex(t2):
                        intsec2.add_unique_point(r)
                elif isinstance(r, geom3d_rat.Segment):
                    if not r.is_triangle_side(t1):
                        intsec1.add_unique_segment(r)
                    if not r.is_triangle_side(t2):
                        intsec2.add_unique_segment(r)
                else:
                    raise Exception('Mesh.delete_self_intersections_rat : complex intersection, '
                                    'not implemented')

            print(f'DSI.Phase.2 : intsec : {i + 1} / {n}')

        if is_log:
            print('DSI.Phase.2 : intsec : end')

        #
        # Third phase. Triangulation and construct outer mesh.
        #

        if is_log:
            print('DSI.Phase.3 : triang : begin')

        # Triangulate triangles.
        for i in range(n):
            (t, intsec) = tis[i]
            if intsec.is_empty():
                m.add_triangle(z, t)
            else:
                small_triangles = t.triangulate(intsec)
                for st in small_triangles:
                    if geom3d_rat.Vector.dot(t.outer_normal, st.outer_normal) < 0:
                        st.flip_normal()
                m.add_triangles(z, small_triangles)
            print(f'DSI.Phase.3 : triang : {i + 1} / {n}')

        if is_log:
            print('DSI.Phase.3 : triang : end')

        # Convert coordinates back to real.
        # We have to convert coordinates only for result mesh
        # because self mesh now is not valid and can not be used.
        m.convert_coordinates_rat_to_real()

        #
        # Fourth phase. Walk outer surface of the mesh.
        #

        if is_log:
            print('DST.Phase.4 : walkin : begin')

        # Set all marks as -1.
        for f in m.faces:
            f.mark = -1

        # We start walking through the mesh from some minimal face.
        f = min(m.faces, key=lambda f: min([n.p[0] for n in f.nodes]))

        # Breadth-first traversal.
        stack = [f]
        while len(stack) > 0:
            f = stack.pop(0)
            f.mark = 1
            for e in f.edges:
                faces_count = len(e.faces)
                if faces_count == 1:
                    nf = None
                elif faces_count == 2:
                    nf = f.neighbour(e)
                elif faces_count == 4:
                    mp = 0.5 * (e.nodes[0].p + e.nodes[1].p)
                    fn = f.outer_normal()
                    i = np.argmax([fn @ (f.third_node(e).p - mp) for f in e.faces])
                    nf = e.faces[i]
                else:
                    raise Exception('Mesh.delete_self_intersections_rat : edge with '
                                    f'unexpected number of incident faces ({faces_count}) is found')
                if not nf is None:
                    if nf.mark < 0:
                        stack.append(nf)

        # Delete all faces with glo_ids < 0.
        m.delete_faces(lambda f: f.mark < 0)

        if is_log:
            print('DST.Phase.4 : walkin : end')

        return m

#===================================================================================================

if __name__ == '__main__':
    start = time.time()
    mesh = Mesh('../data/meshes/tetrahedron_double.dat')
    mesh.print(True, True)
    geom3d_rat.print_statistics()
    print(f'total time : {time.time() - start}')

# ==================================================================================================
