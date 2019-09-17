"""
Utility functions and classes for 3D hydrostatic ocean model
"""
from __future__ import absolute_import
from .firedrake import *
import os
import numpy as np
import sys
from pyop2.profiling import timed_region, timed_function, timed_stage  # NOQA
from mpi4py import MPI  # NOQA
import ufl  # NOQA
import coffee.base as ast  # NOQA
from collections import OrderedDict, namedtuple  # NOQA
from firedrake import Function as FiredrakeFunction
from firedrake import Constant as FiredrakeConstant
from firedrake import Expression as FiredrakeExpression
from abc import ABCMeta, abstractmethod

class FrozenClass(object):
    """A class where creating a new attribute will raise an exception if _isfrozen == True"""
    _isfrozen = False

    def __setattr__(self, key, value):
        if self._isfrozen and not hasattr(self, key):
            raise TypeError('Adding new attribute "{:}" to {:} class is forbidden'.format(key, self.__class__.__name__))
        super(FrozenClass, self).__setattr__(key, value)


class SumFunction(object):
    """
    Helper class to keep track of sum of Coefficients.
    """
    def __init__(self):
        """
        Initialize empty sum.

        get operation returns Constant(0)
        """
        self.coeff_list = []

    def add(self, coeff):
        """
        Adds a coefficient to self
        """
        if coeff is None:
            return
        self.coeff_list.append(coeff)

    def get_sum(self):
        """
        Returns a sum of all added Coefficients
        """
        if len(self.coeff_list) == 0:
            return None
        return sum(self.coeff_list)


class AttrDict(dict):
    """
    Dictionary that provides both self['key'] and self.key access to members.

    http://stackoverflow.com/questions/4984647/accessing-dict-keys-like-an-attribute-in-python
    """
    def __init__(self, *args, **kwargs):
        if sys.version_info < (2, 7, 4):
            raise Exception('AttrDict requires python >= 2.7.4 to avoid memory leaks')
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


class FieldDict(AttrDict):
    """
    AttrDict that checks that all added fields have proper meta data.

    Values can be either Function or Constant objects.
    """
    def _check_inputs(self, key, value):
        if key != '__dict__':
            from firedrake.functionspaceimpl import MixedFunctionSpace, WithGeometry
            if not isinstance(value, (FiredrakeFunction, FiredrakeConstant)):
                raise TypeError('Value must be a Function or Constant object')
            fs = value.function_space()
            is_mixed = (isinstance(fs, MixedFunctionSpace)
                        or (isinstance(fs, WithGeometry)
                            and isinstance(fs.topological, MixedFunctionSpace)))
            if not is_mixed and key not in field_metadata:
                msg = 'Trying to add a field "{:}" that has no metadata. ' \
                      'Add field_metadata entry to field_defs.py'.format(key)
                raise Exception(msg)

    def _set_functionname(self, key, value):
        """Set function.name to key to ensure consistent naming"""
        if isinstance(value, FiredrakeFunction):
            value.rename(name=key)

    def __setitem__(self, key, value):
        self._check_inputs(key, value)
        self._set_functionname(key, value)
        super(FieldDict, self).__setitem__(key, value)

    def __setattr__(self, key, value):
        self._check_inputs(key, value)
        self._set_functionname(key, value)
        super(FieldDict, self).__setattr__(key, value)


ElementContinuity = namedtuple("ElementContinuity", ["horizontal", "vertical"])
"""
A named tuple describing the continuity of an element in the horizontal/vertical direction.

The field value is one of "cg", "hdiv", or "dg".
"""


def element_continuity(ufl_element):
    """Return an :class:`ElementContinuity` instance with the
    continuity of a given element.

    :arg ufl_element: The UFL element to determine the continuity
        of.
    :returns: A new :class:`ElementContinuity` instance.
    """
    elem = ufl_element
    elem_types = {
        'Discontinuous Lagrange': 'dg',
        'Lagrange': 'cg',
        'Raviart-Thomas': 'hdiv',
        'Q': 'cg',
        'DQ': 'dg',
    }

    if isinstance(elem, ufl.finiteelement.mixedelement.MixedElement):
        elem = elem.sub_elements()[0]
    if isinstance(elem, ufl.finiteelement.mixedelement.VectorElement):
        elem = elem.sub_elements()[0]  # take the elem of first component
    if isinstance(elem, ufl.finiteelement.tensorproductelement.TensorProductElement):
        a, b = elem.sub_elements()
        horiz_type = elem_types[a.family()]
        vert_type = elem_types[b.family()]
    elif isinstance(elem, ufl.finiteelement.hdivcurl.HDivElement):
        horiz_type = 'hdiv'
        vert_type = 'hdiv'
    else:
        horiz_type = elem_types[elem.family()]
        vert_type = horiz_type
    return ElementContinuity(horiz_type, vert_type)


def create_directory(path, comm=COMM_WORLD):
    """
    Create a directory on disk

    Raises IOError if a file with the same name already exists.
    """
    if comm.rank == 0:
        if os.path.exists(path):
            if not os.path.isdir(path):
                raise IOError('file with same name exists', path)
        else:
            os.makedirs(path)
    comm.barrier()
    return path


class VelocityMagnitudeSolver(object):
    """
    Computes magnitude of (u[0],u[1],w) and stores it in solution
    """
    def __init__(self, solution, u=None, w=None, min_val=1e-6,
                 solver_parameters={}):
        """
        :arg solution: scalar field for velocity magnitude scalar :class:`Function`
        :type solution: :class:`Function`
        :kwarg u: horizontal velocity
        :type u: :class:`Function`
        :kwarg w: vertical velocity
        :type w: :class:`Function`
        :kwarg float min_val: minimum value of magnitude. Minimum value of solution
            will be clipped to this value
        :kwarg dict solver_parameters: PETSc solver options


        If ``u`` is None computes magnitude of (0,0,w).

        If ``w`` is None computes magnitude of (u[0],u[1],0).
        """
        self.solution = solution
        self.min_val = min_val
        function_space = solution.function_space()
        test = TestFunction(function_space)
        tri = TrialFunction(function_space)

        a = test*tri*dx
        s = 0
        if u is not None:
            s += u[0]**2 + u[1]**2
        if w is not None:
            s += w**2
        l = test*sqrt(s)*dx
        self.prob = LinearVariationalProblem(a, l, solution)
        self.solver = LinearVariationalSolver(self.prob, solver_parameters=solver_parameters)

    def solve(self):
        """Compute the magnitude"""
        self.solver.solve()
        np.maximum(self.solution.dat.data, self.min_val, self.solution.dat.data)



def beta_plane_coriolis_params(latitude):
    r"""
    Computes beta plane parameters :math:`f_0,\beta` based on latitude

    :arg float latitude: latitude in degrees
    :return: f_0, beta
    :rtype: float
    """
    omega = 7.2921150e-5  # rad/s Earth rotation rate
    r = 6371.e3  # Earth radius
    # Coriolis parameter f = 2 Omega sin(alpha)
    # Beta plane approximation f_beta = f_0 + Beta y
    # f_0 = 2 Omega sin(alpha_0)
    # Beta = df/dy|_{alpha=alpha_0}
    #      = (df/dalpha*dalpha/dy)_{alpha=alpha_0}
    #      = 2 Omega cos(alpha_0) /R
    alpha_0 = 2*np.pi*latitude/360.0
    f_0 = 2*omega*np.sin(alpha_0)
    beta = 2*omega*np.cos(alpha_0)/r
    return f_0, beta


def beta_plane_coriolis_function(latitude, out_function, y_offset=0.0):
    """
    Interpolates beta plane Coriolis function to a field

    :arg float latitude: latitude in degrees
    :arg out_function: :class:`Function` where to interpolate
    :kwarg float y_offset: offset (y - y_0) used in Beta-plane approximation.
        A constant in mesh coordinates.
    """
    # NOTE assumes that mesh y coordinate spans [-L_y, L_y]
    f0, beta = beta_plane_coriolis_params(latitude)
    coords = SpatialCoordinate(out_function.function_space().mesh())
    out_function.interpolate(f0 + beta * (coords[1] - y_offset))


def select_and_move_detectors(mesh, detector_locations, detector_names=None,
                              maximum_distance=0.):
    """Select those detectors that are within the domain and/or move them to
    the nearest cell centre within the domain

    :arg mesh: Defines the domain in which detectors are to be located
    :arg detector_locations: List of x, y locations
    :arg detector_names: List of detector names (optional). If provided, a list
       of selected locations and a list of selected detector names are returned,
       otherwise only a list of selected locations is returned
    :arg maximum_distance: Detectors whose initial locations is outside the domain,
      but for which the nearest cell centre is within the specified distance, are
      moved to this location. By default a maximum distance of 0.0 is used, i.e
      no detectors are moved.
    """
    # auxilary function to test whether we can interpolate it in the given locations
    V = FunctionSpace(mesh, "CG", 1)
    v = Function(V)

    P0 = FunctionSpace(mesh, "DG", 0)
    VP0 = VectorFunctionSpace(mesh, "DG", 0)
    dist = Function(P0)
    loc_const = Constant(detector_locations[0])
    xy = SpatialCoordinate(mesh)
    p0xy = Function(VP0).interpolate(xy)

    # comparison operator that sorts on first entry first, etc.
    def min_lexsort(x, y, datatype):
        for xi, yi in zip(x, y):
            if xi < yi:
                return x
            elif yi < xi:
                return y
        # all entries the same:
        return x
    min_lexsort_op = MPI.Op.Create(min_lexsort, commute=False)

    def move_to_nearest_cell_center(location):
        loc_const.assign(location)
        dist.interpolate(dot(xy-loc_const, xy-loc_const))
        ind = dist.dat.data_ro.argmin()
        # smallest distance to a cell centre location on this process:
        local_loc = list(p0xy.dat.data_ro[ind])
        local_dist = np.sqrt(dist.dat.data_ro[ind])
        # select the smallest distance on all processes. If some distances are equal, pick a unique loc. based on lexsort
        global_dist_loc = mesh.comm.allreduce([local_dist]+local_loc, op=min_lexsort_op)
        return global_dist_loc[0], global_dist_loc[1:]

    accepted_locations = []
    accepted_names = []
    if detector_names is None:
        names = [None] * len(detector_locations)
    else:
        names = detector_names
    for location, name in zip(detector_locations, names):
        try:
            v(location)
        except PointNotInDomainError:
            moved_dist, location = move_to_nearest_cell_center(location)
            if moved_dist > maximum_distance:
                continue
        accepted_locations.append(location)
        accepted_names.append(name)

    min_lexsort_op.Free()

    if detector_names is None:
        return accepted_locations
    else:
        return accepted_locations, accepted_names
