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
# Normal recalculation routines (recalc_norms) (for when ctrl-n fails)
# ----------------------------------------------------------------------------

# from mesh_heal.op_norms import *
# recalc_norms(C.active_object)

# Initialization
from .op_gen import *

# ----------------------------------------------------------------------------

class MeshHealRecalcNormsOperator(bpy.types.Operator):
    """Recalculate normals (Mesh Heal)"""
    bl_idname = "mesh.mesh_heal_recalc_norms"
    bl_label = "MH Recalculate Normals"

    @classmethod
    def poll(cls, context):
        ob = context.active_object
        return (ob and ob.type == 'MESH' and context.mode == 'OBJECT')

    def execute(self, context):
        n = recalc_norms(context.active_object)
        self.report({'INFO'}, "%d problem faces selected" % n)
        return {'FINISHED'}


def recalc_norms(obj):
    """Main routine for recalculating normals of object obj
    using a recursive outward face normal casting and propagation method.
    Object is assumed to consist of any number of closed, manifold,
    and potentially intersecting surfaces. Shared edges and vertices 
    are allowed.
    
    Algorithm supports layers of surfaces inside surfaces. 
    Each surface layer normals are set top oppose direction of 
    outer layer normals.

    This method relies on getting correct results from casting 
    phase. Therefore a few heuristics (search for overlapping and 
    high aspect ration faces) are included to choose
    which faces are not included in casting.

    Routine does not guarantee that all face normals are corrected.
    Faces whose normals could not be calculated are selected at the
    end, so the user may choose to operate on those afterwards.
    Overlapping faces are most problematic. Works best for 
    triangles and untwisted and unskewed quadrilateral faces.
    """
    
    # Initialization    
    # Create bmesh
    bm = bmesh_copy_from_object(obj, transform=True, triangulate=False)
    bm.faces.ensure_lookup_table()
    bm.edges.ensure_lookup_table()    
    
    # Face process tracking list
    f_is_processed = len(bm.faces) * [False]

    # Bad face (unsuitable for casting) list
    f_is_noncastable = len(bm.faces) * [False]      

    # First find overlapping faces. Such faces can't be used for casting,
    # since ray casting from such faces can produce wrong results.
    recalc_norms_find_overlapping_neighbor_faces(bm, f_is_noncastable)

    # Add high aspect ratio faces to list of noncastables
    recalc_norms_find_high_aspect_ratio_faces(bm, f_is_noncastable)    
    
    # Deduce normals iteratively until all good faces are processed
    n_unprocessed = f_is_processed.count(False)
    iter = 0
    while n_unprocessed > 0:
        iter += 1
        # First is casting phase
        recalc_norms_cast(obj, bm, f_is_processed, f_is_noncastable)
        
        # Stop if casting does nothing
        if n_unprocessed == f_is_processed.count(False):
            break
        
        # Second is propagation phase
        recalc_norms_propagate(bm, f_is_processed)
        n_unprocessed = f_is_processed.count(False)
        l.debug("Finished main iter %d" % iter \
        + ", unprocessed = %d" % n_unprocessed)

    # Select unprocessed problematic faces
    if f_is_processed.count(False) > 0:
        for i, val in enumerate(f_is_processed):
            if f_is_processed[i]:
                continue
            bm.faces[i].select = True
            l.debug("Selected error face %d" % i)
            
        l.info("Uncomplete finish. Added %d " % n_unprocessed 
        + "unprocessed faces to selection.")
    else:
        l.debug("Complete finish, all normals recalculated!")

    # Save final bmesh back to object and clean up
    bmesh_to_object(obj, bm)
    bm.free()
    del f_is_processed
    return n_unprocessed

def recalc_norms_find_overlapping_neighbor_faces(bm, f_is_overlapping):
    """Finds overlapping neighbor faces. This is done by going through 
    all manifold edges and checking if angle between two connected faces 
    is very small. If yes, face indices are marked as True.
    """

    # MIN_COS_ANGLE is minimum value for cos(angle) below which angle is
    # assumed to be effectively zero
    MIN_COS_ANGLE = 1e-2
    
    for e in bm.edges:
        # Calculate cosine of angle between faces
        cos_epsilon = face_face_cos_angle(e)
        if cos_epsilon == None:
            continue
        elif (cos_epsilon > (1 - MIN_COS_ANGLE)):
            i1 = e.link_faces[0].index
            i2 = e.link_faces[1].index
            f_is_overlapping[i1] = True
            f_is_overlapping[i2] = True
            l.debug("Marked overlapping faces %d and %d" % (i1, i2))

def recalc_norms_find_high_aspect_ratio_faces(bm, f_is_high_aspect_ratio):
    """Finds faces with high aspect ratio. Face aspect ratio is 
    length of longest edge divided by length of shortest
    edge. If a value is exceeded, face index is marked as 
    a high aspect ratio face.
    """

    # AR is maximum allowed face aspect ratio
    AR = 1/0.6 # 60% (or more) of max length is seen as OK
    
    for f in bm.faces:
        min=float_info.max
        max=0.0
        for e in f.edges:
            length = e.calc_length()
            if length < min:
                min = length
            if length > max:
                max = length
        aspect_ratio = max/min
        if aspect_ratio > AR:
            f_is_high_aspect_ratio[f.index] = True
            l.debug("Marked bad aspect ratio face %d" % f.index)
    
def recalc_norms_cast(obj, bm, f_is_processed, f_is_noncastable):
    """Cast phase of recalc_norms: Cast a ray from 
    face center towards normal and opposite direction. 
    If possible, use normal information from what is hit to
    set a normal of this face. Casting is not done for faces
    marked as True in index list f_is_noncastable.
    """
    
    # working copy of processed faces list, to prevent recursive processing
    f_is_processed_work = list(f_is_processed)
    
    n = 0 # number of faces for which normal cast was successful
    
    for i, val in enumerate(f_is_processed):
        # Skip processed faces
        if val == True:
            continue

        # Skip casting for bad faces
        if f_is_noncastable[i]:
            continue

        # Get face f, it's normal and center
        f = bm.faces[i]
        f.normal_update()
        f_no = f.normal
        f_center = f.calc_center_median()

        # Get normal and face index number of what is hit by rays in
        # normal and opposite to normal directions.
        dir_normal, i_normal = \
            obj_probe_hit_face(obj, bm, f_center, f_no, i)
        dir_opposite, i_opposite = \
            obj_probe_hit_face(obj, bm, f_center, -1 * f_no, i)

        l.debug("Cast results for face %d: " % i \
                + "%s %s " % (dir_normal, i_normal) \
                + "%s %s" % (dir_opposite, i_opposite))
        
        # Case 1: both hit world boundary -> do not process
        if (i_normal == -1 and i_opposite == -1):
            continue

        # Case 2: only normal direction hits world boundary -> mark as processed
        elif (i_normal == -1):
            f_is_processed_work[i] = True

        # Case 3: only opposite direction hits world boundary -> flip face normal
        elif (i_opposite == -1):
            f.normal_flip()
            f_is_processed_work[i] = True
        
        # Case 4: both hit unprocessed faces -> do not process        
        elif not f_is_processed[i_normal] and not f_is_processed[i_opposite]:
            continue

        # Case 5: both hit processed faces, process if face normal directions
        # are consinstent
        elif f_is_processed[i_normal] and f_is_processed[i_opposite]:
            if dir_normal!=dir_opposite: 
                f_is_processed_work[i] = True

        # Case 6: face in normal direction has been processed, take normal
        # from there
        elif f_is_processed[i_normal]:
            if dir_normal:
                f_is_processed_work[i] = True
            else:
                f.normal_flip()
                f_is_processed_work[i] = True

        # Case 7: face in opposite direction has been processed, take normal
        # from there
        elif f_is_processed[i_opposite]:
            if dir_opposite:
                f.normal_flip()
                f_is_processed_work[i] = True
            else:
                f_is_processed_work[i] = True
        n += 1
        
    # Write processed states back to original list
    for i, val in enumerate(f_is_processed):
        f_is_processed[i] = f_is_processed_work[i]
        
    l.debug("Cast %d new normals " % n
    + "(total faces %d)" % len(f_is_processed))
            
def recalc_norms_propagate(bm, f_is_processed):
    """Propagation phase of recalc_norms: Copy normals to neighboring faces.
    This is done only for unprocessed faces if at least one neighboring face
    has been processed (so it contains a correct normal). Propagation is 
    repeated until unprocessed face count is not decreasing any more.
    """

    # Set up working copies of processed faces list, to prevent
    # recursive processing
    f_is_processed_work = list(f_is_processed)
    f_is_processed_work2 = list(f_is_processed_work)    
    n_unprocessed = len(f_is_processed)
    iter = 0
    while n_unprocessed > f_is_processed_work.count(False):
        iter += 1
        n_unprocessed = f_is_processed_work.count(False)
        for i, val in enumerate(f_is_processed_work):
            # Skip processed faces
            if val == True:
                continue
            # Try to propagate normal from neighbor to this face
            f_is_processed_work2[i] = \
                recalc_norms_propagate_face(bm, f_is_processed_work, i)
        # Iteration round finished, write intermediate status to
        # work face list
        for i, val in enumerate(f_is_processed_work):
            f_is_processed_work[i] = f_is_processed_work2[i]
        l.debug("Propagation iter %d" % iter \
        + ", unprocessed = %d" % f_is_processed_work.count(False))
        
    # Write processed states back to original list
    for i, val in enumerate(f_is_processed):
        f_is_processed[i] = f_is_processed_work[i]
        
def recalc_norms_propagate_face(bm, f_is_processed, i):
    """Propagates normal to face i from it's neighbor face, if possible. 
    If a neighboring face has been processed, align normal direction
    from it and and return True to mark face i as processed. 
    If no neighbor faces have been processed, return False.
    """
    
    # MIN_COS_ANGLE is clipping value for cos(angle)
    MIN_COS_ANGLE = 1e-2
    
    f = bm.faces[i]
    for e in f.edges:
        l.debug("Processing edge %d" % e.index) 

        # Only manifold edges are considered
        if len(e.link_faces) != 2:
            continue

        # Find neighboring face (f_neighbor)
        if e.link_faces[0] == f:
            f_neighbor = e.link_faces[1]
        else:
            f_neighbor = e.link_faces[0]
        
        # If neighbor has not been processed, go to next edge
        if not f_is_processed[f_neighbor.index]:
            continue
        
        l.debug("Face %d found neighbor face %d" % (i, f_neighbor.index))

        if propagate_face_normal_from_neighbor(f, f_neighbor, e) == True:
            return True
                
    return False

def propagate_face_normal_from_neighbor(f, f_neighbor, e=None):
    """Propagates normal to face f from it's neighbor face f_neighbor,
    connected by edge e. e is searched internally if not given as argument. 
    Returns True if normal propagation was successful or False in all other
    cases (faces are not connected or faces overlap).
    """
    
    # MIN_COS_ANGLE is clipping value for cos(angle)
    MIN_COS_ANGLE = 1e-2

    # Search for edge e or check validity of given edge e
    if e == None:
        e = [e for e in f.edges if e in f_neighbor.edges]
        if e == []:
            return False
        e = e[0]
    elif (e not in f.edges) and (e not in f_neighbor.edges):
        raise ValueError("edge %d does not connect " % e.index \
            + "faces %d and %d" % (f.index, f_neighbor.index))
    
    l.debug("Propagation for face %d " % f.index \
        + "from neighbor face %d " % f_neighbor.index \
        + "using edge %d" % e.index)
    
    # Calculate cos(epsilon), cosine of edge angle between faces
    cos_epsilon = face_face_cos_angle(e, f, f_neighbor)
    if cos_epsilon == None:
        raise ValueError("Calculation of cos_epsilon failed")
                
    # Calculate cos(eta) where eta is the angle between
    # this face current normal and f_neighbor current normal
    cos_eta = f.normal @ f_neighbor.normal

    # Heuristic to flip (or not to flip) this face normal:
    #
    # Case 1: For effectively orthogonal faces (90 degree
    # angle between faces): 
    # calculate cos(theta) where theta is the angle between 
    # f_neighbor_center -> f_center and f.normal, and
    # calculate cos(omega) where omega is the angle between
    # f_neighbor_center -> f_center and f_neighbor.normal, and
    # use those angle cosines to decide normal direction
    if (abs(cos_epsilon) < MIN_COS_ANGLE):
        vec_f2f = f.calc_center_median() - f_neighbor.calc_center_median() 
        vec_f2f.normalize()
        cos_theta = vec_f2f @ f.normal
        cos_omega = vec_f2f @ f_neighbor.normal
        if (cos_theta * cos_omega > 0):
            f.normal_flip()
        l.debug("Orthogonal Face %d, " % f_neighbor.index \
                + "cos_theta=%.8f, " % cos_theta \
                + "cos_omega=%.8f" % cos_omega)
    
    # Case 2: For effectively planar faces (180 degree
    # angle between faces): normals should be almost same
    elif (cos_epsilon < (-1 + MIN_COS_ANGLE)):
        if (cos_eta < 0):
            f.normal_flip()
        l.debug("Planar Face %d, " % f_neighbor.index \
                + "cos_epsilon=%.8f, " % cos_epsilon \
                + "cos_eta=%.8f" % cos_eta)
    
    # Case 3: Do nothing for extremely sharp angles (overlapping faces)
    elif (cos_epsilon > (1 - MIN_COS_ANGLE)):
        l.debug("Overlapping Face %d, " % f_neighbor.index \
                + "cos_epsilon=%.8f, " % cos_epsilon \
                + "cos_eta=%.8f" % cos_eta)
        return False
    
    # Case 4: For other angles, flip this face normal if it is pointing
    # towards wrong direction
    elif (cos_epsilon * cos_eta > 0):
        f.normal_flip()
        l.debug("Flipped Face %d, " % f_neighbor.index \
                + "cos_epsilon=%.8f, " % cos_epsilon \
                + "cos_eta=%.8f" % cos_eta)
    else:
        l.debug("Normal Face %d, " % f_neighbor.index \
                + "cos_epsilon=%.8f, " % cos_epsilon \
                + "cos_eta=%.8f" % cos_eta)
    return True

def propagate_face_normal_from_any(f, avoid_faces=None):
    """Tries to propagate face normal to face f from any one of it's 
    neighboring faces which are not overlapping with face f.
    Do not take normal from faces listed in avoid_faces.
    Returns True if face propagation was succesful and False otherwise.
    """

    for e in f.edges:
        l.debug("Processing edge %d" % e.index) 

        # Only manifold edges are considered
        if len(e.link_faces) != 2:
            continue

        # Find neighboring face (f_neighbor)
        if e.link_faces[0] == f:
            f_neighbor = e.link_faces[1]
        else:
            f_neighbor = e.link_faces[0]
        
        l.debug("Face %d found neighbor face %d" % (f.index, f_neighbor.index))

        # Skip this face if it is listed in avoid_faces
        if avoid_faces and f_neighbor in avoid_faces:
            continue
        
        if propagate_face_normal_from_neighbor(f, f_neighbor, e) == True:
            return True                

    return False
    
