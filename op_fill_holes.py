# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8-80 compliant>

# ----------------------------------------------------------------------------
# Algorithm to find and fill boundary loops in mesh
# ----------------------------------------------------------------------------

# from mesh_heal.op_fill_holes import *
# fill_holes_sharp(C.active_object)

# Initialization
from .op_gen import *
from .op_norms import *

from copy import copy

# ----------------------------------------------------------------------------

class EdgePath:
    """Container class for topology (connectivity information) of edge paths.
    Stores information about path edges, path end points and branches at 
    end points.
    """
    
    def __init__(self, x):
        """Initialize instance with instance x (a BMEdge or another 
        EdgePath)
        """
        
        # Initialization with an edge
        if isinstance(x, bmesh.types.BMEdge):
            self.edges = [] # list of edges in edge path
            # edge, vertex and branch edges at end of edge path
            self.end1_e = x
            self.end1_v = x.verts[0]
            self.end1_branches = [x]
            self.end2_e = x
            self.end2_v = x.verts[1]
            self.end2_branches = [x]
            
        # Initialization with another EdgePath
        elif isinstance(x, EdgePath):
            self.edges = copy(x.edges)
            self.end1_e = x.end1_e
            self.end1_v = x.end1_v
            self.end1_branches = copy(x.end1_branches)
            self.end2_e = x.end2_e
            self.end2_v = x.end2_v
            self.end2_branches = copy(x.end2_branches)

        else:
            raise TypeError("Unknown type")

    def __str__(self):
        str = "EdgePath %s contains\n" % hex(id(self))
        str += " %d edges: " % self.len()
        for i in self.edges:
            str += "%d " % i.index
        str += "\n"

        str += " end1 edge %d," % self.end1_e.index \
            + " end1 vertex %d" % self.end1_v.index \
            + " branches: "
        for i in self.end1_branches:
            str += "%d " % i.index
        str += "\n"
            
        str += " end2 edge %d," % self.end2_e.index \
            + " end2 vertex %d" % self.end2_v.index \
            + " branches: "
        for i in self.end2_branches:
            str += "%d " % i.index
        str += "\n"
        return str
        
    def len(self):
        """Returns number of edges in this edge path"""
        return len(self.edges)

    def extend(self, e, edges):
        """Extends EdgePath towards branch edge e, using edges as
        search domain.
        """


        # Probe at end1
        new_elist = []
        new_vlist = []
        new_branch_edges = []
        
        if e in self.end1_branches:
            n_branches = FHS_get_boundary_edge_path(e, self.end1_v, edges, \
                self.edges, new_elist, new_vlist, new_branch_edges)
            l.debug("new_branch_edges %d" % len(new_branch_edges))
            l.debug("new_branch_edges %s" % hex(id(new_branch_edges)))
            for ne in new_elist:
                if ne not in self.edges:
                    self.edges.append(ne)                    
            self.end1_branches = copy(new_branch_edges)
            self.end1_v = new_vlist[-1]
            if len(new_elist) > 0:
                self.end1_e = new_elist[-1]
            else:
                self.end1_e = []

        # Probe at end2
        new_elist = []
        new_vlist = []
        new_branch_edges = []

        if e in self.end2_branches:
            n_branches = FHS_get_boundary_edge_path(e, self.end2_v, edges, \
                self.edges, new_elist, new_vlist, new_branch_edges)
            l.debug("new_branch_edges %d" % len(new_branch_edges))
            l.debug("new_branch_edges %s" % hex(id(new_branch_edges)))
            for ne in new_elist:
                if ne not in self.edges:
                    self.edges.append(ne)                    
            self.end2_branches = copy(new_branch_edges)
            self.end2_v = new_vlist[-1]
            if len(new_elist) > 0:
                self.end2_e = new_elist[-1]
            else:
                self.end2_e = []

    
                
        l.debug("branch1 %d" % len(self.end1_branches))
        l.debug("branch2 %d" % len(self.end2_branches))
                
        # if (e not in self.end1_branches) and (e not in self.end2_branches):
        #     raise ValueError("Edge %d is not a continuation branch" % e.index)
                
    def is_closed(self):
        """Returns True if edges in this EdgePath form a closed loop,
        otherwise returns False
        """
        return self.end1_v == self.end2_v

# ----------------------------------------------------------------------------

class MeshHealFillHolesSharpOperator(bpy.types.Operator):
    """Fill Holes Sharp (Mesh Heal)"""
    bl_idname = "mesh.mesh_heal_fill_holes_sharp"
    bl_label = "MH Fill Holes (Sharp)"

    @classmethod
    def poll(cls, context):
        ob = context.active_object
        return (ob and ob.type == 'MESH' and context.mode in {'OBJECT','EDIT_MESH'})

    def execute(self, context):
        saved_mode = context.active_object.mode
        n = fill_holes_sharp(context.active_object)
        bpy.ops.object.mode_set(mode = saved_mode)
        self.report({'INFO'}, "%d problem edges selected" % n)
        return {'FINISHED'}

    
def fill_holes_sharp(obj):
    """Fills boundary edges in object obj with triangles by a
    'sharpest angle first' approach. The method processes each 
    continuous boundary edge loop and fills it with
    triangles by connecting edge vertices. Connections are done in 
    sharpness order, so that two edges connected to the vertex with 
    sharpest edge angle are connected by a face. Connections are not 
    made for vertex pairs which are mutually occluded (a neighboring face 
    of the edge loop exists between vertices). Triangle filling is repeated 
    until all boundary edges that can be processed, have been processed.
    """
    
    # Initialization, create bmesh
    bpy.ops.object.mode_set(mode = 'OBJECT')

    # Load original mesh
    bm_orig = bmesh_from_object(obj)
    bm_orig.faces.ensure_lookup_table()
    bm_orig.edges.ensure_lookup_table()
    bm_orig.verts.ensure_lookup_table()

    bm_filled = bmesh.new() # Empty bmesh for saving filled holes
    searched_edges = [] # List of edges that have been searched
    max_edges = len(bm_orig.edges) # maximum number of edges in original mesh

    # Get all boundary edges to be processed
    work_edges = [e for e in bm_orig.edges if len(e.link_faces) == 1]
    
    iter = 0
    while True:
        iter += 1
        l.debug("=====")
        l.info("Loop fill iteration %d, " % iter \
               + "%d work edges remaining" % len(work_edges))

        # Stop if there are no more edges to be processed
        if not work_edges:
            break
        
        # Find boundary edge loop. If edges do not form a closed loop,
        # remove unbranching edges from work edges and go to next edges
        edge_loop = []
        if not FHS_get_continuous_boundary_edges(work_edges, edge_loop):
            work_edges = [e for e in work_edges if e not in edge_loop]
            continue

        # Remove edge loop edges from list of searchable edges
        work_edges = [e for e in work_edges if e not in edge_loop]

        # Stop if edges contains a newly created edge
        # Not needed any more?
        #n_new_edges = [e for e in edges if e.index > max_edges]
        #if len(n_new_edges) > 0:
        #    l.info("No more original boundary edges found! Stopping.")
        #    break

        l.debug("Found %d connected boundary edges" % len(edge_loop))

        # Create bmesh subset of mesh around boundary edges
        bm_sub, edges_sub = FHS_create_subset_mesh(edge_loop)

        # DEBUG: Save subset mesh
        #for e in edges_sub:
        #    e.select = True
        #bmesh_to_object(obj, bm_sub)
        #return 0

        # Create BVHTree of subset mesh for ray casting
        tree = mathutils.bvhtree.BVHTree.FromBMesh(bm_sub)
    
        # Generate vertex data and try to fill holes
        bdata_v, bdata_e1, bdata_e2, bdata_a, bdata_processed, bdata_conn = \
            FHS_generate_vertex_data(tree, edges_sub)
        rval = FHS_fill_holes_to_boundary_edge_path(bm_sub, bm_filled, bdata_v, \
            bdata_e1, bdata_e2, bdata_a, bdata_processed, bdata_conn)

        # Report failed edges
        for e in edges_sub:
            if len(e.link_faces) == 1:
                l.debug("Failed subset edge %d" % e.index)
        l.debug("Failed %d out of %d edges" % (rval, len(edges_sub)))

        del [bdata_v, bdata_e1, bdata_e2, bdata_a, bdata_processed]
        del [bdata_conn, edge_loop]
        bm_sub.free()
        del tree

    # Copy all new faces to original bmesh
    l.info("Merging %d new faces to original mesh.." % len(bm_filled.faces))
    bmesh_copy_face(bm_filled.faces, bm_orig)

    # Select and count remaining boundary edges
    n_unprocessed = 0
    for e in bm_orig.edges:
        if len(e.link_faces) == 1:
            e.select = True
            n_unprocessed +=1
        
    # Save final bmesh back to object and clean up
    bmesh_to_object(obj, bm_orig)
    bm_orig.free()
    bm_filled.free()
    bpy.ops.object.mode_set(mode = 'OBJECT')
    return n_unprocessed

def FHS_create_subset_mesh(edges):
    """Returns a subset bmesh (and subset edges) that includes only 
    faces connected to argument edges.
    """

    bm_sub = bmesh.new()
    
    src_faces = [] # source face list
    vmap = {} # mapping dictionary from old to new vertices
    
    # Find source faces that connect to vertices of edges
    for e in edges:
        for v in e.verts:
            for f in v.link_faces:
                if f in src_faces:
                    continue
                src_faces.append(f)

    # Construct subset mesh consisting of source faces
    for f in src_faces:
        nverts = [] # vertex list for a face
        for v in f.verts:
            if v not in vmap:
                nv = bm_sub.verts.new(v.co)
                bm_sub.verts.index_update()
                bm_sub.verts.ensure_lookup_table()
                # vmap[v] = nv
                vmap[v] = bm_sub.verts[-1]
            else:
                nv = vmap[v]
            nverts.append(nv)
        bm_sub.faces.new(nverts)

    # Update subset mesh lookup tables to prevent index "-1" messing indexing
    bm_sub.faces.index_update()
    bm_sub.edges.index_update()
    bm_sub.verts.index_update()
    bm_sub.faces.ensure_lookup_table()
    bm_sub.edges.ensure_lookup_table()
    bm_sub.verts.ensure_lookup_table()
    
    # Find subset edges that correspond to argument edges
    edges_sub = []
    for es in edges:
        v1 = vmap[es.verts[0]]
        v2 = vmap[es.verts[1]]
        e = [e for e in bm_sub.edges if (v1 in e.verts and v2 in e.verts)][0]
        edges_sub.append(e)

    return bm_sub, edges_sub


def FHS_get_continuous_boundary_edges(work_edges, edge_loop):
    """Finds shortest possible closed boundary edge loop consisting of
    boundary edges, and contains first edge in work_edges. Loop edges are
    appended to list edge_loop. 
    At minimum the first edge is added to edge_loop.
    Return True if closed edge loop was found, and False otherwise.
    """

    if not work_edges:
        raise ValueError("sanity check")
    
    # If there are no branches, give results
    e0 = work_edges[0]
    ep0 = EdgePath(e0)
    ep0.extend(e0, work_edges)
    l.debug(ep0)
    if ep0.is_closed():
        edge_loop += ep0.edges
        del ep0
        return True

    # Otherwise, find combinations of all branches
    eplist = []
    for e in ep0.end1_branches + ep0.end2_branches:
        enew = EdgePath(ep0) # Start with ep0
        enew.extend(e, work_edges) # Extend towards branch
        eplist.append(enew)

    #l.debug("eplist:")
    #for ep in eplist:
    #    l.debug(ep)
        
    # Find first among shortest edge path among closed edge paths:
    loop_found = False
    ep = [ep for ep in eplist if ep.is_closed()]
    if ep:
        loop_found = True
        shortest = min(ep, key=lambda x: x.len())
        edge_loop += shortest.edges

        debugtext = "Chosen loop is "
        for e in edge_loop:
            debugtext += "%d " % e.index
        l.debug(debugtext)

    # if no closed loops were found, return the first unbranched part
    else:
        edge_loop += ep0.edges
        l.debug("No closed loops found for edge %d" % e0.index)
        
    # Cleanup
    for i in ep:
        del i
    del ep0

    return loop_found
    
def FHS_get_boundary_edge_path(e0, v0, edges, old_elist, new_elist, new_vlist, new_branch_edges):
    """Finds edges and vertices that are located along an unbranching edge
    path which starts from vertex v0 towards edge e0.
    Boundary edges are saved to list new_elist and vertices to list new_vlist. 
    Searching is stopped if found edges exist in list old_elist.
    Possible branch edges at the end point are appended to new_branch_edges.
    Only edges in argument list "edges" are eligible for searching.
    Returns number of branches at the end point.
    """

    new_elist.append(e0)
    nv = e0.other_vert(v0) # next vertex
    new_vlist.append(nv)
    
    eold = e0 # old edge
    continue_iteration = True
    
    while continue_iteration:
        # Populate list of boundary edges branching from nv
        test_branches = []
        for e in nv.link_edges:
            if e == eold:
                l.debug("ignore old link_edge %d" %e.index)
                continue
            if (e in old_elist) or (e == e0):
                continue_iteration = False
                l.debug("found old link_edge %d" %e.index)
            if e in edges:
                l.debug("append new approved link_edge %d" %e.index)
                test_branches.append(e)
            else:
                l.debug("not allowed link_edge %d" %e.index)
                
                
        n_branches = len(test_branches)
        l.debug("n_branches %d" % n_branches)
        
        # If there is only one continuation edge, add it and continue iteration
        if continue_iteration and n_branches == 1:
            ne = test_branches[0] # update next edge
            new_elist.append(ne)
            eold = ne
            nv = ne.other_vert(nv) # update next vertex
            new_vlist.append(nv)

        # If branch point or end point was found, stop here
        if n_branches != 1:
            new_branch_edges += test_branches
            l.debug("length of new_branch_edges is %d" % len(new_branch_edges))
            break

    return n_branches


def NOTUSED_FHS_edges_form_closed_loop(bm, edges):
    """Returns True if edges form a closed edge loop in bmesh bm,
    otherwise returns False.
    """

    # dictionary for storing vertex occurrence count
    vdict = {}

    for e in edges:
        for v in e.verts:
            if v not in vdict:
                vdict[v] = 1
            else:
                vdict[v] +=1

    for v in vdict:
        if vdict[v] != 2:
            l.debug("Edge loop is NOT closed")
            return False

    l.debug("Edge loop is closed")
    return True
                
def FHS_generate_vertex_data(source, edges):
    """Generates allowed connectivity data and boundary vertex data
    for vertices of boundary edge path edges for object obj.
    source is a BVHTree or object for ray casting.
    
    Returned boundary vertex data include:
    bdata_v - vertex object 
    bdata_e1 - edge 1 object
    bdata_e2 - edge 2 object
    bdata_a - angle between edge1 and edge2
    bdata_processed - list to mark which vertices have been processed already
    bdata_conn - Numpy matrix of allowed/prohibited connections between 
    vertex pairs. Matrix element is True if ray casting from one vertex 
    to another does not hit a surface in between vertex points, meaning
    that it should be safe to create edge between vertices.    
    """

    bdata_v = []
    bdata_e1 = []
    bdata_e2 = []
    bdata_a = []

    # Go through edges and initialize boundary vertex data
    for e in edges:
        for v in e.verts:
            FHS_populate_vert_data(v, e, bdata_v, bdata_e1, bdata_e2, bdata_a)
    bdata_processed = len(bdata_v) * [False]
    
    # Set up bdata_conn
    nv = len(bdata_v)
    bdata_conn = numpy.full((nv, nv), True, dtype=bool)

    for i in range(nv):
        for j in range(i + 1, nv):
            l.debug("Ray casting: vertex %d to %d" \
                % (bdata_v[i].index, bdata_v[j].index))
            co = bdata_v[i].co # ray start point
            ray_dir = bdata_v[j].co - bdata_v[i].co # ray direction vector
            fi_list = obj_index_list(bdata_v[i].link_faces) # faces of vert i
            ok, hit_co, hit_no, hit_index = \
                face_ray_cast(source, co, ray_dir, fi_list, ray_dir.length)

            # If ray casting hit a surface, then mark connection i,j as 
            # forbidden
            if ok:
                l.debug("Mark verts %d and %d as unconnectable" \
                    % (bdata_v[i].index, bdata_v[j].index))
                bdata_conn[i,j] = False
                bdata_conn[j,i] = False

    return bdata_v, bdata_e1, bdata_e2, bdata_a, bdata_processed, bdata_conn

def FHS_populate_vert_data(v, e, bdata_v, bdata_e1, bdata_e2, bdata_a):
    """Populates vertex data for boundary vertex v which belongs to edge e.
    See FHS_generate_vertex_data for description of other arguments.
    """

    # If vertex is already in boundary data, then add edge and calculate angle
    if v in bdata_v:
        l.debug("Modify data for vertex %d" % v.index)
        ind = bdata_v.index(v)
        bdata_e2[ind] = e
        bdata_a[ind] = \
            calc_vert_edge_edge_angle(v, bdata_e1[ind], bdata_e2[ind])

    # Otherwise, initialize new boundary data for vertex v
    else:
        l.debug("Initialize data for vertex %d" % v.index)
        bdata_v.append(v)
        bdata_e1.append(e)
        bdata_e2.append(None)
        bdata_a.append(math.radians(180))

def FHS_fill_holes_to_boundary_edge_path(bm, bm_filled, bdata_v, \
        bdata_e1, bdata_e2, bdata_a, bdata_processed, bdata_conn):
    """Fills triangle faces to boundary edge path using argument data.
    bm is bmesh used for determining what to fill. Fill faces are
    copied to bm_filled.
    Other arguments are described in function FHS_generate_vertex_data().
    Returns number of failed faces (failed fills).
    """

    # Number of triangle faces to be filled to a boundary loop
    # is number of vertices - 2
    n_unprocessed = bdata_processed.count(False) - 2
    iter = 0
    while n_unprocessed > 0:
        iter += 1
        l.debug("FHS boundary iteration %d " % iter \
            + "unprocessed verts %d" % n_unprocessed)
        # boolean to mark a succesfful face creation
        addition = False
        # Process vertices in increasing angle order
        angles = numpy.argsort(bdata_a)
        for i in angles:
            # Skip if already processed
            if bdata_processed[i] == True:
                continue
            # This vertex and it's edges
            v = bdata_v[i]
            e1 = bdata_e1[i]
            e2 = bdata_e2[i]
            if e1 == None:
                raise ValueError("e1 is None")
            if e2 == None:
                raise ValueError("e2 is None")
                
            l.debug("Processing vertex %d, " % v.index \
                    + "edges %d and %d" % (e1.index, e2.index))
            
            # Get other vertices, if there are edges
            ov1 = None
            if e1 != None:
                ov1 = e1.other_vert(v)
            ov2 = None
            if e2 != None:
                ov2 = e2.other_vert(v)
            # Get next vertex after other vertices
            oov1 = FHS_other_vert(ov1, e1, bdata_v, bdata_e1, bdata_e2)
            oov2 = FHS_other_vert(ov2, e2, bdata_v, bdata_e1, bdata_e2)
            if oov1 == None:
                raise ValueError("oov1 is None")
            if oov2 == None:
                raise ValueError("oov2 is None")
            l.debug("other vertices: %d and %d, " % (ov1.index, ov2.index) \
                + "oo %d and %d" % (oov1.index, oov2.index))

            # Connection between ov1 and ov2 must be allowed
            test1 = FHS_connection_is_allowed(bdata_v, bdata_conn, ov1, ov2)

            # At least one fill direction must remain possible:
            # If a face would be created among vertices v, ov1 and ov2,
            # then ov1 would later on need to be connected with oov2, or
            # ov2 with oov1.
            test2 = FHS_connection_is_allowed(bdata_v, bdata_conn, ov1, oov2)
            test3 = FHS_connection_is_allowed(bdata_v, bdata_conn, ov2, oov1)
            
            # Additionally, connection is not made if it creates sharp edges.
            # This can happen e.g. when boundary has inward extension with
            # three edges connected to the inward vertex.            
            test4 = not FHS_connection_creates_sharp_edge(e1, ov2)
            test5 = not FHS_connection_creates_sharp_edge(e2, ov1)
            
            if not test1:
                l.debug("test1 failed")
            if not test2:
                l.debug("test2 failed")
            if not test3:
                l.debug("test3 failed")
            if not test4:
                l.debug("test4 failed")
            if not test5:
                l.debug("test5 failed")
            
            if test1 and (test2 or test3) and test4 and test5:
                addition = FHS_connect(bm, bm_filled, i, bdata_v, bdata_e1, \
                    bdata_e2, bdata_a)
                bdata_processed[i] = True
                if addition:
                    n_unprocessed -= 1
                break # get out of for loop
        # If no new connections (edges) were made, then stop
        if not addition:
            break

    return n_unprocessed

def FHS_other_vert(v, e, bdata_v, bdata_e1, bdata_e2):
    """Returns the 'other' vertex connected to vertex v than the one on the
    other end of edge e. v is assumed to be connected to boundary edges
    e and e2 in the boundary data. This function finds the vertex at the 
    other end of e2. Returns None in case of failure.
    """

    if v == None or e == None:
        raise ValueError("v or e is None")
    i = bdata_v.index(v)
    l.debug("index %d" % i)
    if bdata_e1[i] == e:
        oe = bdata_e2[i]
        l.debug("oe is %d" % oe.index)
    elif bdata_e2[i] == e:
        oe = bdata_e1[i]
        l.debug("oe is %d" % oe.index)
    else:
        raise ValueError("Can't find other edge")
    if oe == None:
        raise ValueError("oe is None")
    l.debug("other vertices of oe: %d and %d" % (oe.verts[0].index, oe.verts[1].index))
    l.debug("vertex index is %d" % v.index)
    ov = oe.other_vert(v)
    l.debug("found other vertex %d" % ov.index)
    return ov
        
def FHS_connection_is_allowed(bdata_v, bdata_conn, v1, v2):
    """Checks if connection of vertices v1 and v2 is allowed according
    to data in bdata_conn. Vertex index bdata_v is needed for index lookups.
    Return True for allow connection and False for prohibit connection.
    """
    
    # If v1 or v2 are undefined, allow them to be connected.
    # Undefined vertices originate from searching connection beyond
    # loop end points when boundary loop is not closed.
    # TODO: Fix FHS_connect to make this unneeded.
    if v1 == None or v2 == None:
        return True
    
    if v1 in bdata_v:
        iv1 = bdata_v.index(v1)
    else:
        raise ValueError("Vertice %d is not in bdata_v" % v1.index)         

    if v2 in bdata_v:
        iv2 = bdata_v.index(v2)
    else:
        raise ValueError("Vertice %d is not in bdata_v" % v2.index)         
    
    return bdata_conn[iv1, iv2]

def FHS_connection_creates_sharp_edge(e, v):
    """Checks if edge e would create a sharp edge if a
    triangle face were created connecting edge e vertices and vertex v.
    Returns True if edge would become sharp and False otherwise.
    """

    # MIN_COS_ANGLE is minimum value for cos(angle) below which angle is
    # assumed to be effectively zero
    MIN_COS_ANGLE = 1e-2

    if e == None or v == None:
        raise ValueError("e or v is None")

    for f in e.link_faces:
        cos_angle = edge_face_vertex_cos_angle (e, f, v)
        if (cos_angle > (1 - MIN_COS_ANGLE)):
            return True

    return False
            
def FHS_connect(bm, bm_filled, i, bdata_v, bdata_e1, bdata_e2, bdata_a):
    """Creates face to boundary edges connected to vertex with index i
    and updates boundary edge and angle data accordingly.
    Return True if face was created (or exists already), False otherwise.
    """

    v = bdata_v[i]
    e1old = bdata_e1[i]
    e2old = bdata_e2[i]

    # Do nothing for unclosed edge loop end points
    if e1old == None:
        return False
    if e2old == None:
        return False
    
    v1 = e1old.other_vert(v)
    v2 = e2old.other_vert(v)
    
    l.debug("FHS connect starting: ind %d=vert %d, " % (i, v.index) \
        + "edges %d and %d, " % (e1old.index, e2old.index) \
        + "other verts %d and %d" % (v1.index, v2.index))
    
    # Stop if there is already a face that uses these vertices
    flist = bmesh_face_get_partial(bm, [v, v1, v2])
    if flist != None:
        l.debug("Face %d already uses these vertices, " % flist[0].index \
            + "no new face created")
        return True

    # Create face
    neogeo = bm.faces.new([v, v1, v2])
    # Must update face normal manually to make it coherent
    neogeo.normal_update()
    # Must update indices, othewise index is -1 for all new items
    # and that messes up this algorithm.
    bm.edges.index_update()
    bm.faces.index_update()
    bm.faces.ensure_lookup_table()
    bm.edges.ensure_lookup_table()
    # Set face normal direction from any neighbor face.
    # Not necessary here, so commented out.
    # propagate_face_normal_from_any(neogeo)

    # Copy face to fill bmesh
    l.debug(bm_filled)
    bmesh_copy_face(neogeo, bm_filled)
    
    l.debug("Created face %d " % bm.faces[-1].index \
            + "connecting vertices %d " % v.index \
            + "+ %d and %d" % (v1.index, v2.index))

    # Find the edge between v1 and v2. It can be a new edge or a
    # previously existing edge if mesh is non-manifold.
    e = [e for e in bm.edges if v1 in e.verts and v2 in e.verts]
    if not e:
        raise ValueError("sanity check")
    e = e[0]
    l.debug("Last edge index is %d" % e.index)

    # Update new neigbor edge info
    iv1 = bdata_v.index(v1)
    if bdata_e1[iv1] == e1old:
        bdata_e1[iv1] = e
        l.debug("e1 v1 edge update: %d -> %d" % (e1old.index, e.index))
    elif bdata_e2[iv1] == e1old:
        bdata_e2[iv1] = e
        l.debug("e2 v1 edge update: %d -> %d" % (e1old.index, e.index))
    else:
        raise ValueError("Did not find edge %d " % e.index
            + "for v1 %d" % v1.index)
            
    iv2 = bdata_v.index(v2)
    if bdata_e1[iv2] == e2old:
        bdata_e1[iv2] = e
        l.debug("e1 v2 edge update: %d -> %d" % (e2old.index, e.index))
    elif bdata_e2[iv2] == e2old:
        bdata_e2[iv2] = e
        l.debug("e2 v2 edge update: %d -> %d" % (e2old.index, e.index))
    else:
        raise ValueError("Did not find edge %d " % e.index
            + "for v2 %d" % v2.index)
    
    # Recalculate angles. If angle calculation fails: keep calm,
    # make it 180 degrees and carry on.
    bdata_a[iv1] = \
        calc_vert_edge_edge_angle(v1, bdata_e1[iv1], bdata_e2[iv1])
    if bdata_a[iv1] == None:
        bdata_a[iv1] = math.radians(180)

    bdata_a[iv2] = \
        calc_vert_edge_edge_angle(v2, bdata_e1[iv2], bdata_e2[iv2])
    if bdata_a[iv2] == None:
        bdata_a[iv2] = math.radians(180)
    
    return True
