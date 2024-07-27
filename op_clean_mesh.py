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
# Mesh clean up routines (clean_mesh)
# ----------------------------------------------------------------------------

# from mesh_heal.op_clean_mesh import *
# clean_and_patch(C.active_object)

# Initialization
from .op_gen import *
import bmesh

# ----------------------------------------------------------------------------

class MeshHealCleanAndPatchOperator(bpy.types.Operator):
    """Clean and Patch (Mesh Heal). Removes non-manifold surfaces and fills open boundaries"""
    bl_idname = "mesh.mesh_heal_clean_and_patch"
    bl_label = "MH Clean and Patch"

    @classmethod
    def poll(cls, context):
        ob = context.active_object
        return (ob and ob.type == 'MESH' and context.mode == 'OBJECT')

    def execute(self, context):
        saved_mode = context.active_object.mode
        n = clean_and_patch(context.active_object)
        bpy.ops.object.mode_set(mode = saved_mode)
        self.report({'INFO'}, "%d problem verts selected" % n)
        return {'FINISHED'}

class MeshHealSimpleCleanOperator(bpy.types.Operator):
    """Simple Clean (Mesh Heal). Merges closeby vertices, deletes overlapping faces, removes non-manifold elements"""
    bl_idname = "mesh.mesh_heal_simple_clean"
    bl_label = "MH Simple Clean"

    @classmethod
    def poll(cls, context):
        ob = context.active_object
        return (ob and ob.type == 'MESH' and \
            context.mode in {'OBJECT','EDIT_MESH'})

    def execute(self, context):
        saved_mode = context.active_object.mode
        mdist = bpy.context.scene.mesh_heal.vert_merge_distance
        nv0 = len(context.active_object.data.vertices)
        clean_mesh_simple_clean(context.active_object, mdist)
        nv = len(context.active_object.data.vertices)
        bpy.ops.object.mode_set(mode = saved_mode)
        self.report({'INFO'}, "%d verts deleted" % (nv0 - nv))
        return {'FINISHED'}
    

def clean_and_patch(obj):
    """Main mesh cleaning routine. This routine attempts 
    (but does not guarantee) to make closed volumes by
    merging vertices, then repeatedly removing bad faces and 
    refilling boundary holes left in the mesh. 
    
    Current implementation gives up after two iterations. 
    Intersecting faces and boundary edges may remain in mesh.
    """

    # First simple clean
    bpy.ops.mesh.mesh_heal_simple_clean()

    n_verts = 1 # initialize number of bad vertices
    i = 0 # clean-up round counter
    max_iter = 2 # maximum iterations

    while n_verts > 0:
        i += 1
        
        # Give up after maxiumum iterations
        if i > max_iter:
            break
        
        # Delete bad faces and remove dangling edges and verts
        bpy.ops.object.mode_set(mode = 'EDIT')
        clean_mesh_select_bad_verts()
        bpy.ops.mesh.delete(type='ONLY_FACE')
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.delete_loose(use_verts=True, use_edges=True)

        # Fill holes with faces, where possible
        clean_mesh_select_bad_verts()        
        # Sometimes fill connects totally unconnected vertices
        # bpy.ops.mesh.fill()
        # bpy.ops.mesh.edge_face_add() # This doesn't work always right
        # bpy.ops.mesh.quads_convert_to_tris(quad_method='BEAUTY', \
        #                                    ngon_method='BEAUTY')
        # bpy.ops.mesh.fill_holes(sides=0) # not working always
        bpy.ops.mesh.mesh_heal_fill_holes_sharp() # own method, slow
        
        # TODO: Dig down to why Blender's face filling routines fail
        
        # Vertex deletion (or edge collapse) are destructive and 
        # (combined with intersecting faces check) can potentially
        # recursively eat away a lot of geometry. Not good.
        # This is never run, left here just to document this attempt.
        if i > max_iter:
            clean_mesh_select_bad_verts()     
            bpy.ops.mesh.delete(type='VERT') # or edge collapse
            # Fill holes with faces, where possible
            clean_mesh_select_bad_verts()        
            bpy.ops.mesh.fill_holes(sides=0) # or another fill method
            bpy.ops.mesh.quads_convert_to_tris(quad_method='BEAUTY', \
                                               ngon_method='BEAUTY')

        # Count remaining bad vertices
        # BMesh method does not get updates to object data, so force
        # updates via Object mode instead
        # bm = bmesh.from_edit_mesh(obj.data)
        # verts = [ v.index for v in bm.verts if v.select ]
        # bm.free()
        clean_mesh_select_bad_verts()
        bpy.ops.object.mode_set(mode = 'OBJECT')
        verts = [v for v in obj.data.vertices if v.select]
        n_verts = len(verts)
        l.debug("Iteration %d: " % i \
                + "bad verts remaining: %d" % n_verts)

    l.info("Stopped at iteration %d" % i \
           + ", bad verts remaining: %d" % n_verts)
    return n_verts

def obj_select_intersecting_faces(obj):
    """Check and select any faces which intersect.
    """
    import mathutils

    if not obj.data.polygons:
        return None

    bm = bmesh_copy_from_object(obj)
    tree = mathutils.bvhtree.BVHTree.FromBMesh(bm, epsilon=0.00001)
    overlap = tree.overlap(tree)
    faces_error = {i for i_pair in overlap for i in i_pair}
    bm.faces.ensure_lookup_table()
    for i in faces_error:
        bm.faces[i].select_set(True)
    bmesh_to_object(obj, bm)

def clean_mesh_select_bad_verts():
    """Selects bad vertices in mesh. Bad vertices include intersecting and 
    multiple overlapping faces, and boundary vertices. 
    
    Note: Non-contiguous verts are not considered, since they result from
    inconsistent normals. Recalculation of normals later on solves problem.

    Note 2: Shared vertices (or edges) are not considered a problem either. 
    If needed you can use "rip" tool to get rid of them later on.
    """
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(type='VERT')
    bpy.ops.mesh.select_all(action='DESELECT')

    # Select intersecting faces
    bpy.ops.mesh.select_mode(type='FACE')
    obj_select_intersecting_faces(bpy.context.active_object)
    bpy.ops.mesh.select_mode(use_extend=True, use_expand=False, type='VERT')

    # Select also all non-manifold boundaries and faces
    bpy.ops.mesh.select_non_manifold(use_wire=True, use_boundary=True, \
    use_multi_face=True, use_non_contiguous=False, use_verts=False) 
    
def clean_mesh_select_non_manifold_verts():
    """Selects non-manifold vertices in mesh. Non-manifold vertices 
    include vertices in intersecting and multiple overlapping faces, 
    wires and single vertices, but not boundary vertices. 
    """

    bpy.ops.mesh.select_all(action='DESELECT')

    # Select intersecting faces
    bpy.ops.mesh.select_mode(type='FACE')
    obj_select_intersecting_faces(bpy.context.active_object)
    bpy.ops.mesh.select_mode(use_extend=True, use_expand=False, type='VERT')

    # Select also all non-manifold boundaries and faces
    bpy.ops.mesh.select_non_manifold(use_wire=True, use_boundary=False, \
    use_multi_face=True, use_non_contiguous=False, use_verts=False) 

    
def clean_mesh_simple_clean(obj, mdist):
    """Merges closeby vertices, then removes non-manifold vertices, 
    edges and faces from object obj. mdist is vertex merge distance.
    """

    # go to vertex select mode and select all
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_mode(use_extend=False, use_expand=False, type='VERT')
    bpy.ops.mesh.select_all(action='SELECT')

    # initial clean-up, merge closeby vertices
    bpy.ops.mesh.remove_doubles(threshold=mdist)

    # delete overlapping neighbor faces (and dangling edges and verts)
    bpy.ops.mesh.mesh_heal_delete_overlap()
    
    # Delete bad faces
    clean_mesh_select_non_manifold_verts()
    bpy.ops.mesh.delete(type='ONLY_FACE')

    # Delete edges that share >2 faces
    bpy.ops.mesh.select_non_manifold(use_wire=False, use_boundary=False, \
    use_multi_face=True, use_non_contiguous=False, use_verts=False) 
    bpy.ops.mesh.delete(type='EDGE')

    # Remove dangling edges and verts
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.delete_loose(use_verts=True, use_edges=True)
    bpy.ops.object.mode_set(mode='OBJECT')

    
class MeshHealTriangulateTwistedFacesOperator(bpy.types.Operator):
    """Triangulate Twisted Faces (Mesh Heal)"""
    bl_idname = "mesh.mesh_heal_triangulate_twisted"
    bl_label = "MH Triangulate Twists"

    @classmethod
    def poll(cls, context):
        ob = context.active_object
        return (ob and ob.type == 'MESH' and context.mode == 'OBJECT')

    def execute(self, context):
        saved_mode = context.active_object.mode
        n = triangulate_twists(context.active_object)
        bpy.ops.object.mode_set(mode = saved_mode)
        self.report({'INFO'}, "%d faces were split" % n)
        return {'FINISHED'}

def triangulate_twists(ob):
    """Triangulate Twisted Faces (Mesh Heal)"""

    import bmesh
    import math
    bpy.ops.object.mode_set(mode = 'EDIT')
    bm = bmesh.from_edit_mesh(ob.data)

    twistfaces = []

    # Maximum absolute cosine of angle between face normal vector and
    # center-to-corner vector
    angle = bpy.context.scene.mesh_heal.max_abs_twist_angle
    max_abs_cos_alpha = math.cos((90.0 - angle) / 90.0 * math.pi / 2.0)

    # Find all twisted faces
    for f in bm.faces:
        norvec = f.normal
        co = f.calc_center_median()
        for v in f.verts:
            vertvec = v.co - co
            vertvec.normalize()
            abs_cos_alpha = abs(vertvec @ norvec)
            l.debug("face %d abs_cos_alpha %f" % (f.index, abs_cos_alpha))
            if abs_cos_alpha > max_abs_cos_alpha:
                if f not in twistfaces:
                    twistfaces.append(f)
                    f.select = True

    # Triangulate twisted faces
    bmesh.ops.triangulate(bm, faces=twistfaces)

    bmesh.update_edit_mesh(ob.data)
    bm.free()
    return len(twistfaces)

