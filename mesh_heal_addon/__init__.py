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

bl_info = {
    "name": "Mesh Healing Tools",
    "author": "Tuomo Keskitalo",
    "blender": (2, 80, 0),
    "location": "3D View > Sidebar",
    "description": "Utilities for healing arbitrary polygon surface meshes",
    "warning": "Experimental",
    "wiki_url": "https://github.com/tkeskita/mesh_heal",
    "tracker_url": "https://github.com/tkeskita/mesh_heal/issues",
    "support": 'COMMUNITY',
    "category": "Mesh"
}

if "bpy" in locals():
    import importlib
    importlib.reload(op_norms)
    importlib.reload(op_clean_mesh)
    importlib.reload(op_fill_holes)
    importlib.reload(op_sew)
    importlib.reload(op_delete_overlap)
    importlib.reload(op_merge_overlapping_edges)
else:
    import bpy
    from . import (
        op_norms,
        op_clean_mesh,
        op_fill_holes,
        op_sew,
        op_delete_overlap,
        op_merge_overlapping_edges,
        )

from .op_gen import *

class MeshHealSettings(bpy.types.PropertyGroup):
    vert_merge_distance: bpy.props.FloatProperty(
        name="Vertex Merge Distance",
        description="Maximum distance for merging closeby vertices",
        default=0.001,
        precision=5,
        min=0.0, max=float_info.max
    )
    sew_ratio_threshold: bpy.props.FloatProperty(
        name="Max Sew Ratio",
        description="Maximum allowed distance ratio for sewing",
        default=0.3,
        precision=5,
        min=0.001, max=1.0
    )
    max_abs_twist_angle: bpy.props.FloatProperty(
        name="Twist Angle",
        description="Maximum allowed angle (in degrees) for twisted faces",
        default=3.0,
        min=0.0, max=90.0
    )
    max_abs_edge_overlap_angle: bpy.props.FloatProperty(
        name="Edge Overlap Angle",
        description="Maximum allowed angle (in degrees) for determining edge overlapping",
        default=1.0,
        min=0.0, max=90.0
    )
    
class MeshHeal_PT_object_mode(bpy.types.Panel):
    """Mesh Heal tool bar panel in object mode"""
    bl_label = "Mesh Heal (MH)"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "Mesh Heal"
    bl_idname = "VIEW3D_PT_MeshHeal_Object"
    bl_context = "objectmode"

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj and obj.type == 'MESH' and context.mode == 'OBJECT'

    def draw(self, context):
        mesh_heal = context.scene.mesh_heal
        layout = self.layout

        col = layout.column()
        rowsub = col.row(align=True)
        rowsub.operator("mesh.mesh_heal_simple_clean", text="Simple Clean")
        rowsub.prop(mesh_heal, "vert_merge_distance", text="Distance")

        row = layout.row()
        row.operator("mesh.mesh_heal_delete_overlap", \
            text="Delete Overlapping Neighbor Faces")

        col = layout.column()
        rowsub = col.row(align=True)        
        rowsub.operator("mesh.mesh_heal_sew", text="Sew Mesh")
        rowsub.prop(mesh_heal, "sew_ratio_threshold", text="Ratio")

        col = layout.column()
        rowsub = col.row(align=True)
        rowsub.operator("mesh.mesh_heal_triangulate_twisted", text="Triangulate Twists")
        rowsub.prop(mesh_heal, "max_abs_twist_angle", text="Angle")

        row = layout.row()
        row.operator("mesh.mesh_heal_clean_and_patch", text="Clean and Patch")
        row = layout.row()
        row.operator("mesh.mesh_heal_fill_holes_sharp", text="Fill Holes (Sharp)")

        col = layout.column()
        rowsub = col.row(align=True)
        rowsub.operator("mesh.mesh_heal_merge_overlapping_edges", text="Merge Overlapping Edges")
        rowsub.prop(mesh_heal, "max_abs_edge_overlap_angle", text="Angle")

        row = layout.row()
        row.operator("mesh.mesh_heal_recalc_norms", text="MH Recalc Norms")
        
class MeshHeal_PT_edit_mode(bpy.types.Panel):
    """Mesh Heal tool bar panel in edit mode"""
    bl_label = "Mesh Heal (MH)"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "Mesh Heal"
    bl_idname = "VIEW3D_PT_MeshHeal_Edit"
    bl_context = "mesh_edit"

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj and obj.type == 'MESH' and context.mode == 'EDIT_MESH'

    def draw(self, context):
        layout = self.layout
        box = layout.box()
        col = box.column(align=True)
        col.label(text="TODO: Make", icon='ERROR')
        col.label(text="available in edit mode")

        
# def mesh_heal_addition_in_edit_mode(self, context):
#     """Example of addition to existing panel"""
#     user_prefs = context.user_preferences
#     addon_prefs = user_prefs.addons[__package__].preferences

#     layout = self.layout
#     box = layout.box()
#     col = box.column(align=True)
#     col.label("Mesh Heal is not")
#     col.label("available in edit mode")
#     # Then add this to register function:
#     bpy.types.MESH_PT_print3d_mesh.append(mesh_heal_addition_in_edit_mode)

# Registration

classes = (
    MeshHeal_PT_object_mode,
    MeshHeal_PT_edit_mode,
    op_norms.MeshHealRecalcNormsOperator,
    op_clean_mesh.MeshHealCleanAndPatchOperator,
    op_clean_mesh.MeshHealSimpleCleanOperator,
    op_clean_mesh.MeshHealTriangulateTwistedFacesOperator,
    op_fill_holes.MeshHealFillHolesSharpOperator,
    op_sew.MeshHealSewOperator,
    op_delete_overlap.MeshHealDeleteOverlapOperator,
    op_merge_overlapping_edges.MergeOverlappingEdgesOperator,
    MeshHealSettings,
)

def menu_func(self, context):
    self.layout.operator(op_norms.MeshHealRecalcNormsOperator.bl_idname)
    self.layout.operator(op_clean_mesh.MeshHealCleanAndPatchOperator.bl_idname)
    self.layout.operator(op_clean_mesh.MeshHealSimpleCleanOperator.bl_idname)
    self.layout.operator(op_clean_mesh.MeshHealTriangulateTwistedFacesOperator.bl_idname)
    self.layout.operator(op_fill_holes.MeshHealFillHolesSharpOperator.bl_idname)
    self.layout.operator(op_sew.MeshHealSewOperator.bl_idname)
    self.layout.operator(op_delete_overlap.MeshHealDeleteOverlapOperator.bl_idname)
    self.layout.operator(op_merge_overlapping_edges.MergeOverlappingEdgesOperator.bl_idname)

def register():
    for cls in classes:
        bpy.utils.register_class(cls)

    bpy.types.Scene.mesh_heal = \
        bpy.props.PointerProperty(type = MeshHealSettings)
    if "VIEW3D_MT_object_cleanup" in dir(bpy.types):
        bpy.types.VIEW3D_MT_object_cleanup.append(menu_func)
    
def unregister():
    if "VIEW3D_MT_object_cleanup" in dir(bpy.types):
        bpy.types.VIEW3D_MT_object_cleanup.remove(menu_func)
    for cls in classes:
        bpy.utils.unregister_class(cls)

    del bpy.types.Scene.mesh_heal
    
