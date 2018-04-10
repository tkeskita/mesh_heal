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
# clean_mesh(C.active_object)

# Initialization
from .op_gen import *


# ----------------------------------------------------------------------------

class MeshHealCleanMeshOperator(bpy.types.Operator):
    """Clean Mesh (Mesh Heal)"""
    bl_idname = "mesh.mesh_heal_clean_mesh"
    bl_label = "MH Clean Mesh"

    @classmethod
    def poll(cls, context):
        ob = context.active_object
        return (ob and ob.type == 'MESH' and context.mode == 'OBJECT')

    def execute(self, context):
        n = clean_mesh(context.active_object)
        self.report({'INFO'}, "%d problem verts selected" % n)
        return {'FINISHED'}


def clean_mesh(obj):
    """Main mesh cleaning routine. This routine attempts 
    (but does not guarantee) to make closed volumes by
    merging vertices, then repeatedly removing bad faces and 
    refilling boundary holes left in the mesh. 
    
    Current implementation gives up after two iterations. 
    Intersecting faces and boundary edges may remain in mesh.
    """

    # go to vertex select mode and select all
    bpy.ops.object.mode_set(mode = 'EDIT')
    bpy.ops.mesh.select_mode(use_extend=False, use_expand=False, type='VERT')
    bpy.ops.mesh.select_all(action = 'SELECT')

    # initial clean-up, merge closeby vertices
    bpy.ops.mesh.remove_doubles(threshold=0.0004) # 0.0001
    # dissolve zero area faces and zero length edges
    bpy.ops.mesh.dissolve_degenerate()
    bpy.ops.mesh.select_all(action = 'DESELECT')

    n_verts = 1 # initialize number of bad vertices
    i = 0 # clean-up round counter

    while n_verts > 0:
        i += 1
        
        # Give up after two iterations
        if i > 2:
            break
        
        # Delete bad faces and remove dangling edges and verts
        clean_mesh_select_bad_verts()
        bpy.ops.mesh.delete(type='ONLY_FACE')
        bpy.ops.mesh.select_all(action = 'DESELECT')
        bpy.ops.mesh.print3d_clean_isolated()

        # Fill holes with faces, where possible
        clean_mesh_select_bad_verts()        
        # Sometimes fill connects totally unconnected vertices
        # bpy.ops.mesh.fill()
        # bpy.ops.mesh.edge_face_add() # This doesn't work either
        # bpy.ops.mesh.quads_convert_to_tris(quad_method='BEAUTY', \
        #                                    ngon_method='BEAUTY')
        # This does not always work correctly
        bpy.ops.mesh.fill_holes(sides=0)
        bpy.ops.mesh.quads_convert_to_tris(quad_method='BEAUTY', \
                                           ngon_method='BEAUTY')
        # TODO: Dig down why face filling routines fail
        
        # Vertex deletion (or edge collapse) are destructive and 
        # (combined with intersecting faces check) can potentially
        # recursively eat away a lot of geometry. Not good.
        # This is left here just to document this attempt.
        if i > 2:
            clean_mesh_select_bad_verts()     
            bpy.ops.mesh.delete(type='VERT') # or edge collapse
            # Fill holes with faces, where possible
            clean_mesh_select_bad_verts()        
            # bpy.ops.mesh.fill()
            # bpy.ops.mesh.edge_face_add()
            bpy.ops.mesh.fill_holes(sides=0)
            bpy.ops.mesh.quads_convert_to_tris(quad_method='BEAUTY', \
                                               ngon_method='BEAUTY')

        # Count remaining bad vertices
        clean_mesh_select_bad_verts()   
        bm = bmesh.from_edit_mesh(obj.data)
        verts = [ v.index for v in bm.verts if v.select ]
        n_verts = len(verts)
        l.debug("Iteration %d: " % i \
                + "bad verts remaining: %d" % n_verts)

    l.info("Stopped at iteration %d" % i \
           + ", bad verts remaining: %d" % n_verts)
    bm.free()
    bpy.ops.mesh.select_all(action = 'DESELECT')
    bpy.ops.object.mode_set(mode = 'OBJECT')
    return n_verts

def clean_mesh_select_bad_verts():
    """Selects bad vertices in mesh. Bad vertices include intersecting and 
    multiple overlapping faces, and boundary vertices. 
    
    Note: Non-contiguous verts are not considered, since they result from
    inconsistent normals. Recalculation of normals later on solves problem.

    Note 2: Shared vertices (or edges) are not considered a problem either. 
    If needed you can use "rip" tool to get rid of them later on.
    """
    bpy.ops.mesh.select_all(action = 'DESELECT')

    # Select intersecting faces
    bpy.ops.mesh.select_mode(type='FACE')
    bpy.ops.mesh.print3d_check_intersect() # Find intersecting faces
    bpy.ops.mesh.print3d_select_report() # Select intersecting faces
    bpy.ops.mesh.select_mode(use_extend=True, use_expand=False, type='VERT')

    # Select also all non-manifold boundaries and faces
    bpy.ops.mesh.select_non_manifold(use_wire=True, use_boundary=True, \
    use_multi_face=True, use_non_contiguous=False, use_verts=False) 