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
# sew_mesh(C.active_object)

# Initialization
from .op_gen import *
import mathutils

# ----------------------------------------------------------------------------

class MeshHealSewOperator(bpy.types.Operator):
    """Sew Mesh (Mesh Heal)"""
    bl_idname = "mesh.mesh_heal_sew"
    bl_label = "MH Sew Mesh"

    @classmethod
    def poll(cls, context):
        ob = context.active_object
        return (ob and ob.type == 'MESH' and context.mode == 'OBJECT')

    def execute(self, context):
        threshold = bpy.context.scene.mesh_heal.sew_ratio_threshold
        n = sew_mesh(context.active_object, threshold)
        if n != None:
            self.report({'INFO'}, "Merged %d vertices" % n)
        else:
            self.report({'INFO'}, "Nothing to sew")
        return {'FINISHED'}

    
def sew_mesh(obj, threshold):
    """Reduce number of boundary edges in mesh by merging pairs of closeby 
    boundary vertices (sew open seams). Boundary vertices are merged 
    if ratio of vertex pair distance to smallest boundary edge length
    (of edges connected to either vertex) is smaller than threshold.
    Return number of merged vertices.

    - If vertices are already connected at boundary path and they are
    two hops apart (one other vertex exists in between), merge must be
    prohibited if it creates sharp face angle at in-between vertex. 
    - Must prohibit merge if ray cast shows a blocking face..?
    """

    # Select non-manifold boundaries
    bpy.ops.object.mode_set(mode='EDIT')    
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.mesh.select_mode(use_extend=False, use_expand=False, type='EDGE')
    bpy.ops.mesh.select_non_manifold(use_wire=False, use_boundary=True, \
    use_multi_face=False, use_non_contiguous=False, use_verts=False) 
    bpy.ops.object.mode_set(mode='OBJECT')

    # Initialization, create bmesh
    bm = bmesh_copy_from_object(obj)
    bm.faces.ensure_lookup_table()
    bm.edges.ensure_lookup_table()
    bm.verts.ensure_lookup_table()

    # Get selected boundary vertices and edges to a list
    vlist = [v for v in bm.verts if v.select == True]
    elist = [e for e in bm.edges if e.select == True]
    nv = len(vlist)
    if nv < 5:
        return None

    l.info("Sew using ratio threshold %f" % threshold)
    l.info("Starting neighbor search for %d boundary vertices.." % nv)
    sew_candidates = get_sew_candidates(bm, vlist, threshold)
    l.info("Processing %d vertex pairs.." % len(sew_candidates))
    for i, pair in enumerate(sew_candidates):
        l.debug("Sew candidate #%d: " % i + str(pair))
    n_merges = merge_sew_candidates(bm, vlist, elist, sew_candidates)
    
    # Save final bmesh back to object and clean up
    bmesh_to_object(obj, bm)
    bm.free()
    bpy.ops.object.mode_set(mode='OBJECT')
    return n_merges

def get_sew_candidates(bm, vlist, threshold):
    """Returns vertex pair list of sew candidates from bmesh bm 
    and boundary vertex list vlist.
    threshold is maximum distance ratio allowed for merging.
    Entry [v1, v2] in returned list means that v1 can be merged to v2. 
    If also entry [v2, v1] exists, then also opposite merge is OK. 
    """

    # Generate KD tree for spatial neighbor search
    nv = len(vlist)
    kd = mathutils.kdtree.KDTree(nv)
    for i, v in enumerate(vlist):
        kd.insert(v.co, i)    
    kd.balance()
    
    sew_candidates = [] # list of vertex pairs eligible to be sewn together
    
    for v in vlist:
        # Initialize
        v_c_dist = float_info.max # distance to unconnected vertex
        c = None # closest unconnected vertex
        l.debug("Processing vertex " + str(v))

        # Find unprocessed closest unconnected vertex
        for (co, index, dist) in kd.find_n(v.co, 4):
            vtest = vlist[index]
            l.debug("  KD found vertex %d, dist %f" % (vtest.index, dist))
            if v == vtest:
                l.debug("    Vertex is source vertex, ignore")
                continue
            if bmesh_verts_share_edge(bm, v, vtest)[0] == False:
                v_c_dist = dist
                c = vtest
                l.debug("    Vertex is unconnected, " \
                        + "set c_dist %f" % dist)
                break

        if c == None:
            l.debug("  No unconnected vertices were found")
            continue

        # Calculate shortest distance from v to connected vertices
        v_dist_min = min([e.calc_length() for e in v.link_edges])
        
        # If length ratio is acceptable, add cadidate vertex pair to list
        # along with ratio information
        ratio = v_c_dist / v_dist_min
        if ratio <= threshold:
            sew_candidates.append([ratio, v, c])
            l.debug("  Connection possible, added " \
                + str(v) + " and " + str(c))
        else:
            l.debug("  Unconnected vertex is too far away, no addition")

    del kd
    return sew_candidates

def merge_sew_candidates(bm, vlist, elist, sew_candidates):
    """Merges vertex pairs in sew_candidates in bmesh bm.
    vlist is boundary vertex list and elist boundary edge list.
    Returns number of merged vertices.
    """

    processed = [] # list of processed vertices
    n_merges = 0 # number of merges done
    
    # Merges are done according to increasing ratio value
    sorted_candidates = sorted(sew_candidates, key = lambda x: x[0])
    
    for ratio, v, c in sorted_candidates:
        if v in processed or c in processed:
            continue

        # Disallow merging of one-hop-apart neighbors, because
        # that would create overlapping face
        if is_one_hop_neighbor(v, c, vlist, elist):
            continue
        
        # Move to mid point if both candidates are listed
        if ([c, v] in [x[1:3] for x in sew_candidates]):
            mergepoint = (v.co + c.co) / 2.0
            v.co = mergepoint
            c.co = mergepoint

        # Otherwise move vertex on top of neighbor
        else:
            v.co = c.co

        # mark both vertices as processed
        n_merges += 1
        processed.append(v)
        processed.append(c)
        
    # Merge overlapping vertices
    DIST = 1e-6 # tolerance for overlapping vertex coordinates
    bmesh.ops.remove_doubles(bm, verts=vlist, dist=DIST)
    
    return n_merges

def is_one_hop_neighbor(v, c, vlist, elist):
    """Returns True if vertice v and c are connected by a neighbor
    vertex in vlist, connected by edges listed in edge list elist.
    """

    evlist = [e for e in elist if v in e.verts]
    eclist = [e for e in elist if c in e.verts]
    # Find if there is an vertex in vlist, whose edges are part of
    # v edge list and c edge list (so there is one-hop connection).
    shared_vert = [v for v in vlist if \
        (v.link_edges[0] in evlist and v.link_edges[1] in eclist) or \
        (v.link_edges[1] in evlist and v.link_edges[0] in eclist)]
                   
    if shared_vert:
        return True
    return False
