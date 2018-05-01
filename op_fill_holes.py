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
# Algorithm to fill boundary loops in meshes.
# ----------------------------------------------------------------------------

# from mesh_heal.op_fill_holes import *
# fill_holes_sharp(C.active_object)

# Initialization
from .op_gen import *
from .op_norms import *

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
        n = fill_holes_sharp(context.active_object)
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
    
    processed_edges = [] # List of processed edges
    max_edges = len(bm_orig.edges) # maximum number of edges in original mesh
    iter = 0
    while True:
        iter += 1
        l.info("=====")
        l.info("Loop fill iteration %d" % iter)
        edges = FHS_get_continuous_boundary_edges(bm_orig, processed_edges)
        if edges == []:
            l.info("No unprocessed boundary edges found! Stopping.")
            break

        # Stop if edges contains a newly created edge
        n_new_edges = [e for e in edges if e.index > max_edges]
        if len(n_new_edges) > 0:
            l.info("No more original boundary edges found! Stopping.")
            break

        # Skip if they do not for a closed loop
        if not FHS_edges_form_closed_loop(bm_orig, edges):
            l.debug("Edges do not form a closed loop, skip.")
            continue
        
        l.debug("Found %d connected boundary edges" % len(edges))

        # Create subset of mesh and boundary edges of subset mesh
        bm_sub, edges_sub = FHS_create_subset_mesh(edges)

        # Create BVHTree of subset mesh for ray casting
        tree = mathutils.bvhtree.BVHTree.FromBMesh(bm_sub)
    
        # Generate vertex data and try to fill holes
        bdata_v, bdata_e1, bdata_e2, bdata_a, bdata_processed, bdata_conn = \
            FHS_generate_vertex_data(tree, edges_sub)
        rval = FHS_fill_holes_to_boundary_edge_path(bm_sub, bm_filled, bdata_v, \
            bdata_e1, bdata_e2, bdata_a, bdata_processed, bdata_conn)
        for e in edges_sub:
            if len(e.link_faces) == 1:
                l.debug("Failed edge %d" % e.index)
        l.info("Failed %d out of %d edges" % (rval, len(edges)))

        del [bdata_v, bdata_e1, bdata_e2, bdata_a, bdata_processed]
        del [bdata_conn, edges]
        bm_sub.free()
        del tree


    # Copy all new faces to original bmesh
    for f in bm_filled.faces:
        bmesh_copy_face(f, bm_orig)

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


def FHS_get_continuous_boundary_edges(bm, processed_edges):
    """Finds first continuous set of boundary edges in bmesh bm
    and returns a list of edge objects that are part of that
    set of boundary edges. Note: Edge list is not ordered.
    """
    
    # start the list of boundary edges
    boundary_edges = []

    for e in bm.edges:
        # only boundary edges connected to one and only one face are considered
        if len(e.link_faces) != 1:
            continue
        
        # do not process same edges twice
        if e in processed_edges:
            continue

        # add boundary edge to lists        
        l.debug("Adding boundary edge %d" % e.index)
        processed_edges.append(e)
        boundary_edges.append(e)
        
        # recursively process edge vertices to populate edge list
        for v in e.verts:
            FHS_get_next_boundary_edge(v, e, boundary_edges, processed_edges)
            
        break
    
    return boundary_edges

def FHS_edges_form_closed_loop(bm, edges):
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
                
def FHS_get_next_boundary_edge(v0, e0, boundary_edges, processed_edges):
    """Finds next boundary edge from edge e0 which is connected
    to edge's vertex v0, adds it to boundary_edges list and calls 
    itself for the next edge, to recursively find next edges
    in the boundary path.
    """

    l.debug("Probe next boundary for v %d" % v0.index + " e %d" % e0.index)
    
    candidate = None # candidate as next boundary edge
    
    for e in v0.link_edges:
        if e == e0:
            continue
        if len(e.link_faces) != 1:
            continue
        # Stop if this is the second new boundary edge at this
        # vertex. Boundary branching is not supported.
        if candidate != None:
            return None
        candidate = e

    # Exit if nothing was found
    if candidate == None:
        return None

    # Exit if this edge has been already added
    if candidate in boundary_edges:
        return None
    
    # At this point a new candidate edge has been found, so        
    # add it to list
    l.debug("Added new boundary edge %d" % candidate.index)
    boundary_edges.append(candidate)
    processed_edges.append(candidate)

    # then find the continuation vertex and recurse
    for v in candidate.verts:
        if v == v0:
            continue
        FHS_get_next_boundary_edge(v, candidate, boundary_edges, \
            processed_edges)
        
    return None

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
    """Populates vertex data for boundary vertice v which belongs to edge e.
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
            l.debug("Processing vertex %d, " % v.index \
                    + "edges %d and %d" % (e1.index, e2.index))
            
            # Get other vertices, if there are edges
            ov1 = None
            if e1 != None:
                ov1 = e1.other_vert(v)
            ov2 = None
            if e2 != None:
                ov2 = e2.other_vert(v)
            # Get next vertice after other vertices
            oov1 = FHS_other_vert(ov1, e1, bdata_v, bdata_e1, bdata_e2)
            oov2 = FHS_other_vert(ov2, e2, bdata_v, bdata_e1, bdata_e2)
            l.debug("other vertices: %d and %d," % (ov1.index, ov2.index) \
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
        return None
    i = bdata_v.index(v)
    if bdata_e1[i] == e:
        oe = bdata_e2[i]
    elif bdata_e2[i] == e:
        oe = bdata_e1[i]
    else:
        return None
    if oe == None:
        return None
    return oe.other_vert(v)
        
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
    """Checks if boundary edge e would create a sharp edge if a
    triangle face were created connecting edge e vertices and vertex v.
    Returns True if edge would become sharp and False otherwise.
    """

    # MIN_COS_ANGLE is minimum value for cos(angle) below which angle is
    # assumed to be effectively zero
    MIN_COS_ANGLE = 1e-2

    if e == None or v == None:
        return None

    # only boundary edges connected to one and only one face are considered
    if len(e.link_faces) != 1:
        raise ValueError("Edge %d contains " % e.index \
            + "%d faces" % len(e.link_faces))
    f = e.link_faces[0]
    cos_angle = edge_face_vertex_cos_angle (e, f, v)
    if cos_angle:
        return (cos_angle > (1 - MIN_COS_ANGLE))
    return None
            
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
    bmesh_copy_face(neogeo, bm_filled)
    
    l.debug("Created face %d " % bm.faces[-1].index \
            + "connecting vertices %d " % v.index \
            + "+ %d and %d" % (v1.index, v2.index))

    # If a new edge was created, it is the last one
    e = bm.edges[-1]
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
