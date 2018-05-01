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
# General and auxiliary functions

# ----------------------------------------------------------------------------
# Initialization

import importlib
import bpy
import bmesh
import mathutils
import array
import numpy
import math
from sys import float_info

# Set up logging of messages using logging
# Logging is nicely explained in:
# https://code.blender.org/2016/05/logging-from-python-code-in-blender/
# Note to self: To see debug messages, configure logging in file
# $HOME/.config/blender/{version}/scripts/startup/setup_logging.py
# add there something like:
# import logging
# logging.basicConfig(format='%(funcName)s: %(message)s', level=logging.DEBUG)
import logging as l

# Uses bmesh_copy_from_object() and bmesh_to_object() from
# Print 3D Tools add-on. Use import instead of exec.
# exec(open(bpy.utils.resource_path('LOCAL') + \
#          "scripts/addons/object_print3d_utils/mesh_helpers.py").read())

from object_print3d_utils.mesh_helpers import *

# ----------------------------------------------------------------------------

def print_face_index():
    """Prints index numbers of currently selected faces."""

    obj = bpy.context.active_object
    bm = bmesh_from_object(obj)
    for f in bm.faces:
        if f.select:
            print(f.index)

def print_edge_index():
    """Prints index numbers of currently selected edges."""

    obj = bpy.context.active_object
    bm = bmesh_from_object(obj)
    for f in bm.edges:
        if f.select:
            print(f.index)
            
def print_vertex_index():
    """Prints index numbers of currently selected faces."""

    obj = bpy.context.active_object
    bm = bmesh_from_object(obj)
    for v in bm.verts:
        if v.select:
            print(v.index)

def select_face_index(ilist):
    """Add argument face index or list of face indices to active selection."""
    
    obj = bpy.context.active_object
    bm = bmesh_from_object(obj)
    bm.edges.ensure_lookup_table()    

    # Convert to list if needed, to make it iterable
    if not isinstance(ilist, list):
        ilist = [ilist]

    for i in ilist:    
        bm.faces[i].select = True
        
    bmesh.update_edit_mesh(obj.data, True)    

def select_edge_index(ilist):
    """Add argument edge index or list of edge indices to active selection."""
    
    obj = bpy.context.active_object
    bm = bmesh_from_object(obj)
    bm.edges.ensure_lookup_table()    

    # Convert to list if needed, to make it iterable
    if not isinstance(ilist, list):
        ilist = [ilist]

    for i in ilist:    
        bm.edges[i].select = True
        
    bmesh.update_edit_mesh(obj.data, True)    

def select_vertex_index(ilist):
    """Add argument vertex index or list of vertex indices to 
    active selection.
    """
    
    obj = bpy.context.active_object
    bm = bmesh_from_object(obj)
    bm.verts.ensure_lookup_table()    

    # Convert to list if needed, to make it iterable
    if not isinstance(ilist, list):
        ilist = [ilist]

    for i in ilist:    
        bm.verts[i].select = True
        
    bmesh.update_edit_mesh(obj.data, True)    
    
def obj_index_list(obj_list):
    """Generates index list from argument object list. Use this to convert 
    list of objects to a list of object indices. Index is assumed to be
    returned by object index property. Works for e.g. BMesh faces, edges and
    vertices.
    """

    r = len(obj_list) * [0]
    for i in range(len(r)):
        r[i] = obj_list[i].index
    return r


def bmesh_edge_center(edge):
    """Calculates coordinates for edge center using edge vertex coordinates."""
    return (edge.verts[0].co + edge.verts[1].co) / 2
        
def bmesh_face_get_partial(bm, vlist):
    """Finds faces which uses all vertices in vertex list vlist,
    but may contain also other vertices. Returns None if none is found,
    otherwise returns list of faces.
    """
    
    f = [f for f in bm.faces if all(v in f.verts for v in vlist)]
    if f == []:
        return None
    return f

def obj_probe_hit_face(obj, bm, f_co, f_dir, f_index):
    """Casts a ray towards f_dir from coordinates f_co and
    checks if the normal of the face hit by ray is towards this face.
    f_index is the face index number of source face at f_co, used for 
    self collision check.
    
    Return value 1: True if normal of hit face is towards this face, 
    otherwise False.

    Return value 2: index number of face which was hit by ray casting.
    """

    ok, hit_co, hit_no, hit_index = face_ray_cast(obj, f_co, f_dir, f_index)

    # If ray hits world boundary, then normal direction is right
    if not ok:
        return (True, hit_index)
    
    # Calculate alignment of hit normal direction with this face normal.
    # Note: Can't use hit_no because object data normals can be out of date.
    cos_theta = f_dir * bm.faces[hit_index].normal    
    if cos_theta > 0:
        return (False, hit_index)
    else:
        return (True, hit_index)

def face_ray_cast(source, co, ray_dir, fi_list, max_ray_len=float_info.max):
    """Face based wrapper for mathutils.bvhtree.ray_cast() (source is 
    a BVHTree) and bpy.types.Object.ray_cast() (source is an Object).  
    Guarantees that ray cast from co towards ray_dir will not hit
    face(s) with index in list of indices fi_list.
    This is used to prevent self-collision for ray casting.
    """

    # Convert fi_list to list if needed, to make it iterable
    if not isinstance(fi_list, list):
        fi_list = [fi_list]
        
    # EPS is length of ray_dir that co is projected as a starting 
    # point used in ray casting.
    EPS0 = 1e-6
    EPS = EPS0
    
    ray_start = co

    # ray_cast returns slighty different things for Objects and BVHTrees
    if isinstance(source, bpy.types.Object):
        ok, hit_co, hit_no, hit_index = \
            source.ray_cast(ray_start + EPS*ray_dir, ray_dir, max_ray_len)
    elif isinstance(source, mathutils.bvhtree.BVHTree):
        hit_co, hit_no, hit_index, hit_length = \
            source.ray_cast(ray_start + EPS*ray_dir, ray_dir, max_ray_len)
        if hit_co:
            ok = True
        else:
            ok = False
    else:
        raise TypeError("source type is not compatible")
            
    # Twisted faces can cause ray hitting itself.
    # If that happens, then cast another ray starting from near hit point.
    iter = 0
    maxiter = 12
    while hit_index in fi_list and iter < maxiter:
        ray_start = hit_co
        EPS *= 2 # Increase EPS to speed up convergence
        if isinstance(source, bpy.types.Object):
            ok, hit_co, hit_no, hit_index = \
                source.ray_cast(ray_start + EPS*ray_dir, ray_dir, max_ray_len)
        elif isinstance(source, mathutils.bvhtree.BVHTree):
            hit_co, hit_no, hit_index, hit_length = \
                source.ray_cast(ray_start + EPS*ray_dir, ray_dir, max_ray_len)
            if hit_co:
                ok = True
            else:
                ok = False
        iter += 1
    
    if iter == maxiter:
        raise ValueError("Cast ray reached maximum iterations %d" % maxiter)

    if ok:
        # If hit coordinates are close to target coordinates,
        # then don't count that as a hit.
        # TODO: Check why BVHTree doesn't seem to honor max ray length?
        if (hit_co - co).length / max_ray_len > (1 - EPS0):
            l.debug("Max ray length limit reached")
            return False, None, None, None
    
    return ok, hit_co, hit_no, hit_index
    

def calc_vert_edge_edge_angle(v, e1, e2):
    """Calculates angle between edges e1 and e2, both of which are
    connected at vertex v. Returns None if edges are not connected
    at vertex v, or if any argument is None.
    """
    
    if v == None or e1 == None or e2 == None:
        return None
    if v not in e1.verts or v not in e2.verts:
        return None

    # Generate normalized edge vectors
    if e1.verts[0] == v:
        vec1 = e1.verts[1].co - v.co
    else:
        vec1 = e1.verts[0].co - v.co
    vec1.normalize()

    if e2.verts[0] == v:
        vec2 = e2.verts[1].co - v.co
    else:
        vec2 = e2.verts[0].co - v.co
    vec2.normalize()
    
    # Calculate angle [rad] between vec1 and vec2
    # Limit -1.0 <= vecprod <= 1.0 so that value is physical
    vecprod = max(-1.0, min(1.0, vec1 * vec2))
    angle = math.acos(vecprod) 
    l.debug("Angle between edges %d and %d " % (e1.index, e2.index) \
        + "is %f" % angle)
    return angle

def face_face_cos_angle(e, f=None, f_neighbor=None):
    """Calculates cosine of angle between two connected faces, which 
    share a common edge e. If faces f1 and f2 are not given as arguments,
    it is assumed that edge e is shared by exactly two faces
    (non-manifold edges are not allowed). Returns None if edge does not 
    connect two faces.
    """

    if f == None or f_neighbor == None:
        # Find faces of edge e
        l.debug("Edge %d has %d link_faces" % (e.index, len(e.link_faces)))
        if len(e.link_faces) != 2:
            return None
        f = e.link_faces[0]
        f_neighbor = e.link_faces[1]
    else:
        # Make sure edge connects given two faces before continuing 
        etest = [etest for etest in f.edges if etest in f_neighbor.edges]
        if e not in etest:
            return None
                    
    # face center coordinates
    f_center = f.calc_center_median()
    f_neighbor_center = f_neighbor.calc_center_median()
    return edge_vec_vec_cos_angle(e, f_center, f_neighbor_center)

def edge_vec_vec_cos_angle(e, vec1, vec2):
    """Calculates cosine of angle between two vectors which are
    orthogonal projections for edge e. This is used to calculate
    cos(angle) for two faces which share edge e, and vec1 and vec2
    represent the faces' center coordinate vectors.
    """
                    
    # Calculate cos(epsilon) where epsilon is the angle between the two 
    # faces connected by edge e. cos(epsilon) is mathematically angle between
    # vectors e_ortho_f -> f_center and e_ortho_f_neighbor -> f_neighbor_center.
    # e_ortho_f is point along edge e which forms 90 degree angle between
    # edge point 1 -> edge point 2 and e_ortho_f -> f_center.
    # Similarly e_ortho_f_neighbor is point along edge e which forms 90 degree
    # angle between edge point 1 -> edge point 2 and 
    # e_ortho_f_neighbor -> f_neighbor_center.
    #
    # Diagram example with two triangular faces: + marks face boundary edges,
    # x = e_ortho_f, y = e_ortho_f_neighbor, c = edge or face center
    #
    #      +++
    #     +   +++++ 
    #    +  c       +++++            <-- f
    #   +   I            +++++++
    #  1----x------c---y--------2    <-- edge e
    #   +++++          I      ++
    #        +++++     c   +++       <-- f_neighbor
    #             ++++   ++
    #                 +++
    
    e_center = bmesh_edge_center(e) # edge center
    vec_e = e.verts[1].co - e.verts[0].co # edge vector

    # Calculate orthogonal vector e_ortho_f -> f_center and normalize it
    f_center = vec1 # face center coordinates
    vec_f_center = f_center - e_center 
    project_f = vec_f_center.project(vec_e) # project vec_f to vec_e
    e_ortho_f = e_center + project_f # coordinates for x
    vec_ortho_f = f_center - e_ortho_f # orthogonal vector 
    vec_ortho_f.normalize() # normalize it
    
    # Similarly to above, calculate orthogonal vector 
    # e_ortho_f_neighbor -> f_neighbor_center and normalize it
    f_neighbor_center = vec2
    vec_f_neighbor_center = f_neighbor_center - e_center 
    project_f_neighbor = vec_f_neighbor_center.project(vec_e)
    e_ortho_f_neighbor = e_center + project_f_neighbor
    vec_ortho_f_neighbor = f_neighbor_center - e_ortho_f_neighbor
    vec_ortho_f_neighbor.normalize()

    # Finally calculate cos(angle) between faces
    cos_epsilon = vec_ortho_f * vec_ortho_f_neighbor
    # Limit -1.0 <= cos_epsilon <= 1.0 for physical correctness
    cos_epsilon = max(-1.0, min(1.0, cos_epsilon))
    return cos_epsilon


def edge_face_vertex_cos_angle (e, f, v):
    """Calculates cos(angle) for angle which would be formed
    between face f (which contains edge e) and a hypothetical connected 
    triangle face f2, which is hypothetically formed by connecting vertices
    of edge e to vertex v. Alternatively, v can be thought of as a vertex
    located at the hypothetical face center. Result is same in both cases.
    """

    if e == None or f == None or v == None:
        return None

    if f not in e.link_faces:
        return None

    # face center coordinates
    f_center = f.calc_center_median()
    return edge_vec_vec_cos_angle(e, f_center, v.co)
    

def bmesh_verts_share_edge(bm, v1, v2):
    """Checks if vertices v1 and v2 in bmesh bm are connected
    via an edge. First return value is boolean which is True if
    edge exists and False otherwise. Second return value is the
    shared edge or None.
    """

    if v1 == v2:
        raise ValueError(str(v1) + " and " + str(v2) + " are equal")
    for e in bm.edges:
        if (v1 in e.verts) and (v2 in e.verts):
            return True, e
    return False, None

def bmesh_copy_face(src_flist, bm):
    """Creates a copy of faces in face list src_flist to bmesh bm."""

    from .op_norms import propagate_face_normal_from_any
    
    vlist = [] # list of face vertices in destination bmesh
    flist = [] # list of new faces in destination bmesh
    
    # Convert flist to list if needed, to make it iterable
    if not isinstance(src_flist, list):
        src_flist = [src_flist]

    for f in src_flist:
        for v in f.verts:
            test, vdst = bmesh_vert_exists_at(bm, v.co)
            if test:
                if vdst in vlist:
                    raise ValueError("Vertex %d alredy in list" % v.index)
                vlist.append(vdst)
            else:
                nv = bm.verts.new(v.co)
                vlist.append(nv)

        nf = bm.faces.new(vlist)
        nf.normal_update()
        flist.append(nf)

    # Finally update lookup tables and propagate normals from neighbors
    bm.faces.index_update()
    bm.faces.ensure_lookup_table()

    for f in flist:
        propagate_face_normal_from_any(nf)

def bmesh_vert_exists_at(bm, co):
    """Finds a vertex in bm located very close to coordinates co.
    First return value is True if vertex is found and False otherwise.
    Second return value is the vertex, or None if none was found.
    """

    # tolerance allowed for difference in vertex location
    EPS = 1e-6
    
    # Decided to use absolute tolerance to test for equality
    v = [v for v in bm.verts if (v.co - co).length < EPS]
    # relative error tolerance alternative would be something like
    # if 2.0 * (v.co - co).length / (v.co.length + co.length) < EPS]
    if len(v) > 0:
        return True, v[0]
    return False, None
