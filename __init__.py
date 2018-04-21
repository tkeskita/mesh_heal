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
    "blender": (2, 79, 0),
    "location": "3D View > Toolbox",
    "description": "Utilities for healing arbitrary surface meshes",
    "warning": "WIP",
    # "wiki_url": "https://none",
    "support": 'COMMUNITY',
    "category": "Mesh"
}

if "bpy" in locals():
    import importlib
    importlib.reload(op_norms)
    importlib.reload(op_clean_mesh)
    importlib.reload(op_fill_holes)
    importlib.reload(op_sew)
else:
    import bpy
    from . import (
        op_norms,
        op_clean_mesh,
        op_fill_holes,
        op_sew,
        )

from .op_gen import *

# requirement: 3D Print Toolbox Add-on. Force enable addon if not enabled.

import addon_utils
addon='object_print3d_utils'
loaded_default, loaded_state = addon_utils.check(addon)
l.debug("loaded_default %s" % loaded_default)
l.debug("loaded_state %s" % loaded_state)
if not loaded_default:
    l.info("Enabling %s addon" % addon)
    addon_utils.enable(addon, default_set=True, persistent=True)

class MeshHealSettings(bpy.types.PropertyGroup):
    sew_ratio_threshold = bpy.props.FloatProperty(
        name="Max Sew Ratio",
        description="Maximum allowed distance ratio for sewing",
        default=0.3,
        precision=5,
        min=0.001, max=1.0
    )
    
class MeshHealToolBarInObjectMode(bpy.types.Panel):
    """Mesh Heal tool bar panel in object mode"""
    bl_label = "Mesh Heal (MH)"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    bl_category = "3D Printing"
    bl_idname = "Mesh_heal_object"
    bl_context = "objectmode"

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj and obj.type == 'MESH' and context.mode == 'OBJECT'

    def draw(self, context):
        mesh_heal = context.scene.mesh_heal
        layout = self.layout
        row = layout.row()
        row.operator("mesh.mesh_heal_clean_mesh", text="MH Clean Mesh")
        row = layout.row()
        row.operator("mesh.mesh_heal_fill_holes_sharp", text="MH Fill Holes (Sharp)")
        row = layout.row()
        row.operator("mesh.mesh_heal_recalc_norms", text="MH Recalc Norms")
        col = layout.column()
        rowsub = col.row(align=True)        
        rowsub.operator("mesh.mesh_heal_sew", text="MH Sew Mesh")
        rowsub.prop(mesh_heal, "sew_ratio_threshold", text="ratio")
        
class MeshHealToolBarInEditMode(bpy.types.Panel):
    """Mesh Heal tool bar panel in edit mode"""
    bl_label = "Mesh Heal (MH)"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    bl_category = "3D Printing"
    bl_idname = "Mesh_heal_edit"
    bl_context = "mesh_edit"

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj and obj.type == 'MESH' and context.mode == 'EDIT_MESH'

    def draw(self, context):
        layout = self.layout
        box = layout.box()
        col = box.column(align=True)
        col.label(text="Mesh Heal is not", icon='ERROR')
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
    MeshHealToolBarInObjectMode,
    MeshHealToolBarInEditMode,
    op_norms.MeshHealRecalcNormsOperator,
    op_clean_mesh.MeshHealCleanMeshOperator,
    op_fill_holes.MeshHealFillHolesSharpOperator,
    op_sew.MeshHealSewOperator,
    MeshHealSettings,
)
    
def register():
    for cls in classes:
        bpy.utils.register_class(cls)

    bpy.types.Scene.mesh_heal = \
        bpy.props.PointerProperty(type = MeshHealSettings)
    
def unregister():
    for cls in classes:
        bpy.utils.unregister_class(cls)

    del bpy.types.Scene.mesh_heal
    
