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
        n = sew_mesh(context.active_object, 0.5)
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
    bpy.ops.object.mode_set(mode = 'EDIT')    
    bpy.ops.mesh.select_all(action = 'DESELECT')
    bpy.ops.mesh.select_mode(use_extend=False, use_expand=False, type='EDGE')
    bpy.ops.mesh.select_non_manifold(use_wire=False, use_boundary=True, \
    use_multi_face=False, use_non_contiguous=False, use_verts=False) 
    bpy.ops.object.mode_set(mode = 'OBJECT')

    # Initialization, create bmesh
    bm = bmesh_copy_from_object(obj)
    bm.faces.ensure_lookup_table()
    bm.edges.ensure_lookup_table()
    bm.verts.ensure_lookup_table()

    # Get selected boundary vertices to a list
    vlist = [v for v in bm.verts if v.select == True]
    nv = len(vlist)
    if nv < 5:
        return None

    l.info("Starting candidate search for %d boundary vertices.." % nv)
    sew_candidates = get_sew_candidates(bm, vlist, threshold)
    l.info("Merging %d vertex pairs.." % len(sew_candidates))
    for i, pair in enumerate(sew_candidates):
        l.debug("Sew candidate #%d: " % i + str(pair))
    n_merges = merge_sew_candidates(bm, vlist, sew_candidates)
    
    # Save final bmesh back to object and clean up
    bmesh_to_object(obj, bm)
    bm.free()
    bpy.ops.object.mode_set(mode = 'OBJECT')
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
    processed_verts = [] # list of vertices in sew_candidates
    
    for v in vlist:
        # Skip if already processed
        if v in processed_verts:
            continue
        
        # Initialize
        v_c_dist = float_info.max # distance to unconnected vertex
        c = None # closest unconnected vertex
        l.debug("Processing vertex " + str(v))

        # 1. Find unprocessed closest unconnected vertex
        for (co, index, dist) in kd.find_n(v.co, 4):
            vtest = vlist[index]
            l.debug("  KD found vertex %d, dist %f" % (vtest.index, dist))
            if v == vtest:
                l.debug("    Vertex is source vertex, ignore")
                continue
            if bmesh_verts_share_edge(bm, v, vtest)[0] == False \
               and vtest not in processed_verts:
                v_c_dist = dist
                c = vtest
                l.debug("    Vertex is unprocessed and unconnected, " \
                        + "set c_dist %f" % dist)
                break

        if c == None:
            l.debug("  No unconnected vertices were found")
            continue

        # Calculate shortest distance from v to connected vertices
        v_dist_min = min([e.calc_length() for e in v.link_edges])
        
        # If length ratio is acceptable, add cadidate vertex pair to list
        # and processed list to make each vertex processed only once
        if (v_c_dist / v_dist_min) <= threshold:
            sew_candidates.append([v, c])
            processed_verts.append(v)
            processed_verts.append(c)            
            l.debug("  Connection possible, added " \
                + str(v) + " and " + str(c))
        else:
            l.debug("  Unconnected vertex is too far away, no addition")

    del kd
    return sew_candidates

def merge_sew_candidates(bm, vlist, sew_candidates):
    """Merges vertex pairs in sew_candidates in bmesh bm.
    Returns number of merged vertices.
    """

    # Move vertices to their merging location
    n_merges = 0
    for v, c in sew_candidates:
        # Move to mid point if both candidates are listed
        if ([c, v] in sew_candidates) and (v.co != c.co):
            mergepoint = (v.co + c.co) / 2.0
            v.co = mergepoint
            c.co = mergepoint
            n_merges += 1
        # Otherwise move to last point if only one vertex is listed
        else:
            v.co = c.co
            n_merges += 1

    # Merge overlapping vertices
    DIST = 1e-6 # tolerance for overlapping vertex coordinates
    bmesh.ops.remove_doubles(bm, verts=vlist, dist=DIST)
    
    return n_merges

# PROBLEMS:
# - first in list steals merge from closest
# - sharp edge issue remains
