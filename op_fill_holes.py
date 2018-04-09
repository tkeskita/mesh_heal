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
# Algorithms to fill boundary holes in meshes. Try these if Blender's
# official fill routines don't work for you.
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
        return (ob and ob.type == 'MESH' and context.mode == 'OBJECT')

    def execute(self, context):
        n = fill_holes_sharp(context.active_object)
        self.report({'INFO'}, "%d problem edges selected" % n)
        return {'FINISHED'}

    
def fill_holes_sharp(obj):
    """Fills boundary edges in object obj with triangles by a
    'sharpest angle first' approach. The method processes each 
    continuous boundary edge loop (open or closed) and fills it with
    triangles by connecting edge vertices. Connections are done in 
    sharpness order, so that two edges connected to the vertex with 
    sharpest edge angle are connected by a face. Connections are not 
    made for vertex pairs which are mutually occluded (face exists 
    between vertices. Filling is repeated until all boundary edges that 
    can be processed, have been processed.
    """
    
    # Initialization, create bmesh
    bpy.ops.object.mode_set(mode = 'OBJECT')
    
    bm = bmesh_copy_from_object(obj) # CHECMKE: bmesh_from_object?
    bm.faces.ensure_lookup_table()
    bm.edges.ensure_lookup_table()
    bm.verts.ensure_lookup_table()

    # List of processed edges
    processed_edges = []
    
    iter = 0
    while True:
        iter += 1
        edges = FHS_get_continuous_boundary_edges(bm, processed_edges)
        if edges == []:
            l.info("Loop fill iteration %d" % iter \
                + ", no unprocessed boundary edges found!")
            break
        l.debug("Boundary %d" % iter \
            + ", found %d connected boundary edges" % len(edges))
        l.debug(edges)
        # if len(edges) == 1:
        #    l.debug("Single edge %d is left unprocessed" % edges[0].index)
        #    continue            

        # Generate vertex data and try to fill holes
        bdata_v, bdata_e1, bdata_e2, bdata_a, bdata_processed, bdata_conn = \
            FHS_generate_vertex_data(obj, edges)
        rval = FHS_fill_holes_to_boundary(bm, bdata_v, bdata_e1, bdata_e2, \
            bdata_a, bdata_processed, bdata_conn)
        l.info("Boundary %d:" % iter \
            + " Failed %d, %d edges" % (rval, len(edges)))        
        del [bdata_v, bdata_e1, bdata_e2, bdata_a, bdata_processed]
        del [bdata_conn, edges]
        
    # Save final bmesh back to object and clean up
    bmesh_to_object(obj, bm)
    bm.free()
    bpy.ops.object.mode_set(mode = 'OBJECT')
    return rval
    
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

def FHS_generate_vertex_data(obj, edges):
    """Generates allowed connectivity data and boundary vertex data
    for vertices of boundary edge path edges for object obj.
    
    Returned boundary vertex data include:
    bdata_v - vertex index number
    bdata_e1 - edge1 index number
    bdata_e2 - edge2 index number
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
    FRAC = 0.999 # fraction of length of ray direction for ray tracing

    for i in range(nv):
        for j in range(i + 1, nv):
            l.debug("Ray casting: vertex %d to %d" \
                % (bdata_v[i].index, bdata_v[j].index))
            f_co = bdata_v[i].co # ray start point
            f_dir = bdata_v[j].co - bdata_v[i].co # ray direction vector
            f_index = obj_index_list(bdata_v[i].link_faces) # faces of vert i
            ok, hit_co, hit_no, hit_index = \
                face_ray_cast(obj, f_co, f_dir, f_index, f_dir.length*FRAC)

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

def FHS_fill_holes_to_boundary(bm, bdata_v, bdata_e1, bdata_e2, \
        bdata_a, bdata_processed, bdata_conn):
    """Fills triangle faces to boundary edge path using argument data.
    Arguments are described in function FHS_generate_vertex_data().
    """

    n_unprocessed = bdata_processed.count(False)
    iter = 0
    n_unallowed = 0 # number of rejected face creation attempts
    while n_unprocessed > 0:
        iter += 1
        l.debug("FHS boundary iteration %d " % iter \
            + "unprocessed edges %d" % n_unprocessed)
        # boolean to mark addition of a new edge when creating face
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

            # If a face would be created between vertices v, ov1 and ov2,
            # then ov1 would later on need to be connected with oov2, and
            # ov2 with oov1. Therefore face creation is allowed only if those
            # connections are allowed according to connectivity matrix.
            # Also connection between ov1 and ov2 must be allowed.
            if FHS_connection_is_allowed(bdata_v, bdata_conn, ov1, oov2) and \
                FHS_connection_is_allowed(bdata_v, bdata_conn, ov2, oov1) and \
                FHS_connection_is_allowed(bdata_v, bdata_conn, ov1, ov2):
                    addition = FHS_connect(bm, i, bdata_v, bdata_e1, \
                        bdata_e2, bdata_a)
                    bdata_processed[i] = True
                    n_unprocessed -= 1
                    break # get out of for loop
            else:
                # Face creation was not allowed
                n_unallowed += 1
                
        # If no new connections (edges) were made, then stop
        if not addition:
            break
    return n_unallowed

def FHS_other_vert(v, e, bdata_v, bdata_e1, bdata_e2):
    """Returns the 'other' vertex connected to vertex v than the one on the
    other end of edge e. v is assumed to be connected to boundary edges
    e and e2 in the boundary data. This function finds the vertex at the 
    other end of e2.
    """

    if v == None or e == None:
        return None
    i = bdata_v.index(v)
    if bdata_e1[i] == e:
        oe = bdata_e2[i]
    elif bdata_e2[i] == e:
        oe = bdata_e1[i]
    else:
        raise ValueError ("Could not find other edge for vertex %d" % i)
    if oe == None:
        return None
    return oe.other_vert(v)
        
def FHS_connection_is_allowed(bdata_v, bdata_conn, v1, v2):
    """Checks if connection of vertices v1 and v2 is allowed according
    to data in bdata_conn. Vertex index bdata_v is needed for index lookups.
    Return True for allow connection and False for prohibit connection.
    """
    
    # If v1 or v2 are undefined, allow their hypothetical connection.
    # Undefined vertices originate frot searching connection beyond
    # loop end points when boundary loop is not closed.
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
    
            
def FHS_connect(bm, i, bdata_v, bdata_e1, bdata_e2, bdata_a):
    """Creates face to boundary edges connected to vertex with index i
    and updates boundary edge and angle data accordingly.
    Returns True if boundary loop is not fully closed, and False otherwise. 
    In other words, False indicates that last triangle in boundary edge was 
    filled.
    """

    v = bdata_v[i]
    e1old = bdata_e1[i]
    e2old = bdata_e2[i]

    # Do nothing for unclosed edge loop end points
    if e1old == None:
        return True
    if e2old == None:
        return True
    
    v1 = e1old.other_vert(v)
    v2 = e2old.other_vert(v)
    
    l.debug("FHS connect starting: ind %d=vert %d, " % (i, v.index) \
        + "edges %d and %d, " % (e1old.index, e2old.index) \
        + "other verts %d and %d" % (v1.index, v2.index))
    
    # If there is already a face that uses these vertices, then stop
    flist = bmesh_face_get_partial(bm, [v, v1, v2])
    if flist != None:
        l.debug("Face %d already uses these vertices, " % flist[0].index \
            + "no new face created")
        return True

    # If there already exists edges among vertices, then mark this as last fill
    lastFill = False
    if bm.edges.get([v, v1]) != None and \
        bm.edges.get([v1, v2]) != None and \
        bm.edges.get([v2, v]) != None:
            lastFill = True
    
    neogeo = bm.faces.new([v, v1, v2])
    # Must update face normal manually to make it coherent
    neogeo.normal_update()
    
    # Must update indices, othewise index is -1 for all new items
    # and that messes up this algorithm. TODO: optimize somehow?
    bm.edges.index_update()
    bm.faces.index_update()
    bm.faces.ensure_lookup_table()
    bm.edges.ensure_lookup_table()

    propagate_face_normal_from_any(neogeo)
    
    l.debug("Created face %d " % bm.faces[-1].index \
            + "connecting vertices %d " % v.index \
            + "+ %d and %d" % (v1.index, v2.index))

    # If a new edge was created, it is the last one
    e = bm.edges[-1]
    l.debug("Last edge index is %d" % e.index)

    # Check if this was the last triangle to be filled. If so, stop.
    if lastFill:
        l.debug("Last triangle in this boundary loop was filled!")
        return False

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
    
    # Recalculate angles
    bdata_a[iv1] = \
        calc_vert_edge_edge_angle(v1, bdata_e1[iv1], bdata_e2[iv1])
    bdata_a[iv2] = \
        calc_vert_edge_edge_angle(v2, bdata_e1[iv2], bdata_e2[iv2])
    
    return True
