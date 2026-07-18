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
# Algorithm to merge overlapping edges
# ----------------------------------------------------------------------------

# Initialization
from .op_gen import *
import mathutils

# ----------------------------------------------------------------------------

class MergeOverlappingEdgesOperator(bpy.types.Operator):
    """Merge Overlapping Edges (Mesh Heal). Finds connected boundary edges forming a sharp angle and merges the edges"""
    bl_idname = "mesh.mesh_heal_merge_overlapping_edges"
    bl_label = "MH Merge Overlapping Edges"

    @classmethod
    def poll(cls, context):
        ob = context.active_object
        return (ob and ob.type == 'MESH' and context.mode == 'OBJECT')

    def execute(self, context):
        n = merge_edges(context.active_object)
        if n != None:
            self.report({'INFO'}, "Merged %d edges" % n)
        else:
            self.report({'INFO'}, "No overlapping edges found in object")
        return {'FINISHED'}


def merge_edges(obj):
    """Merge ovelapping edges in selection. Overlap is defined by
    edge-edge angle being very small for a vertex in the
    selection. Overlap is removed by cutting the longer edge at a
    distance from the shorter edge, and merging the shorted edge
    vertex at the new vertex.
    """

    import math
    # Initialization, create bmesh
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.ops.object.mode_set(mode = 'EDIT')
    bm = bmesh.from_edit_mesh(obj.data)

    # Edge-edge angle (converted from degrees to radians) for considering overlap
    small_angle = bpy.context.scene.mesh_heal.max_abs_edge_overlap_angle
    small_angle *= math.pi / 180.0

    # Do edge merges until no more merges are possible
    merge_done = True
    n_merges = 0
    while merge_done:
        bm.normal_update()
        bm.faces.index_update()
        bm.verts.index_update()
        bm.edges.index_update()
        bm.faces.ensure_lookup_table()
        bm.edges.ensure_lookup_table()
        bm.verts.ensure_lookup_table()

        # Find boundary edges
        edges = [e for e in bm.edges if len(e.link_faces) == 1]

        # Find vertices shared by two boundary edges
        vees = []
        for e1 in edges:
            for v in e1.verts:
                e2s = [e for e in edges if e != e1 and e in v.link_edges]
                for e2 in e2s:
                    vees.append((v, e1, e2))

        # Find first vertex edge edge angle which is small enough to
        # indicate overlapping edge exists
        merge_done = False
        for v, e1, e2 in vees:
            angle = calc_vert_edge_edge_angle(v, e1, e2)
            if angle < small_angle:
                bm = merge_edge(bm, v, e1, e2)
                merge_done = True
                n_merges += 1
                break

    # Save final bmesh back to object and clean up
    bmesh.update_edit_mesh(mesh=obj.data)
    bm.free()
    bpy.ops.object.mode_set(mode='OBJECT')
    return n_merges


def merge_edge(bm, v, e1, e2):
    """Merge shorter edge with longer edge of the two edges connected at
    vertex v.
    """
    from sys import float_info
    e1_length = e1.calc_length()
    e2_length = e2.calc_length()
    # l.debug("Processing vertex %d and edges %d and %d" % (v.index, e1.index, e2.index))

    # Cut longer edge (el) by length from shorter edge
    if e1_length > e2_length:
        el = e1
        es = e2
        fac = e2_length / e1_length
    else:
        el = e2
        es = e1
        fac = e1_length / e2_length
    ne, nv = bmesh.utils.edge_split(el, v, fac)

    # Merge shorter edge end vertex with the new vertex
    l.debug("Merging edge vertices %d and %d" % (es.other_vert(v).index, nv.index) + " located at %r and %r" % (es.other_vert(v).co, nv.co))
    el.select = True
    bmesh.ops.remove_doubles(bm, verts=[es.other_vert(v), nv], dist=float_info.max)
    return bm
