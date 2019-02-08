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
# Mesh clean up routines (clean_overlap)
# ----------------------------------------------------------------------------

# Initialization
from .op_gen import *
import itertools

# ----------------------------------------------------------------------------

class MeshHealDeleteOverlapOperator(bpy.types.Operator):
    """Delete Overlap (Mesh Heal)"""
    bl_idname = "mesh.mesh_heal_delete_overlap"
    bl_label = "MH Delete Overlapping Neighbor Faces"

    @classmethod
    def poll(cls, context):
        ob = context.active_object
        return (ob and ob.type == 'MESH' and context.mode in {'OBJECT','EDIT_MESH'})

    def execute(self, context):
        saved_mode = context.active_object.mode
        bpy.ops.object.mode_set(mode = 'OBJECT')
        n = delete_overlap(context.active_object)
        bpy.ops.object.mode_set(mode = saved_mode)
        self.report({'INFO'}, "%d faces deleted" % n)
        return {'FINISHED'}


def delete_overlap(obj):
    """Deletes overlapping neighbor faces and removes dangling edges 
    and vertices. Neighbor faces mean faces that share an edge.
    Overlapping is determined by face-face-edge angle.
    """

    bm = bmesh_copy_from_object(obj, transform=False, triangulate=False)
    bm.faces.ensure_lookup_table()
    bm.edges.ensure_lookup_table()
    bm.verts.ensure_lookup_table()

    # MIN_COS_ANGLE is minimum value for cos(angle) below which angle is
    # assumed to be effectively zero
    MIN_COS_ANGLE = 1e-2

    flist = [] # list of overlapping faces to be deleted
    
    for e in bm.edges:
        if len(e.link_faces) < 2:
            continue
        
        for i, j in itertools.combinations(e.link_faces, 2):
            l.debug("Checking faces %d and %d" % (i.index, j.index))
            # Calculate cosine of angle between faces
            cos_epsilon = face_face_cos_angle(e, i, j)
            if cos_epsilon == None:
                raise ValueError("cos(epsilon) is None")
            l.debug("cos(epsilon) is %f" % cos_epsilon)
            if (cos_epsilon > (1 - MIN_COS_ANGLE)):
                if i in flist or j in flist:
                    continue
                # Primarily remove face with fewer vertices
                if len(i.verts) < len(j.verts):
                    add_to_flist(i, flist)
                    l.debug("first has less verts")
                elif len(i.verts) > len(j.verts):
                    add_to_flist(j, flist)
                    l.debug("second has less verts")
                # Secondarily remove smaller face
                elif i.calc_area() < j.calc_area():
                    add_to_flist(i, flist)
                    l.debug("first has smaller area")
                else:
                    add_to_flist(j, flist)
                    l.debug("second has smaller area")
                
    n_appends = len(flist)    

    # Delete overlappings
    bmesh.ops.delete(bm, geom=flist, context='FACES_ONLY') # context 3 = DEL_ONLYFACES

    # Save final bmesh back to object
    bmesh_to_object(obj, bm)
    bm.free()

    # Remove dangling edges and vertices
    bpy.ops.mesh.print3d_clean_isolated()

    return n_appends

def add_to_flist(f, flist):
    """Selects and appends face f to face list flist, 
    if face is not there already.
    """
    if f in flist:
        return None
    flist.append(f)
    f.select = True
    l.debug("Appending face %d" % f.index)

