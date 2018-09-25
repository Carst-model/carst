"""
Utility functions and classes
"""
from __future__ import absolute_import
from .firedrake import *
import os
import numpy as np
import sys
from .physical_constants import physical_constants
from pyop2.profiling import timed_region, timed_function, timed_stage  # NOQA
from mpi4py import MPI  # NOQA
import ufl  # NOQA
import coffee.base as ast  # NOQA
from collections import OrderedDict, namedtuple  # NOQA
from .field_defs import field_metadata
from .log import *
from firedrake import Function as FiredrakeFunction
from firedrake import Constant as FiredrakeConstant
from abc import ABCMeta, abstractmethod

ds_surf = ds_t
ds_bottom = ds_b


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
            is_mixed = (isinstance(fs, MixedFunctionSpace) or
                        (isinstance(fs, WithGeometry) and
                         isinstance(fs.topological, MixedFunctionSpace)))
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
    }

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


def comp_volume_2d(eta, bath):
    """Computes volume of the 2D domain as an integral of the elevation field"""
    val = assemble((eta+bath)*dx)
    return val


class VelocityMagnitudeSolver(object):
    """
    Computes magnitude of (u[0],u[1]) and stores it in solution
    """
    def __init__(self, solution, u=None, min_val=1e-6,
                 solver_parameters={}):
        """
        :arg solution: scalar field for velocity magnitude scalar :class:`Function`
        :type solution: :class:`Function`
        :kwarg u: horizontal velocity
        :type u: :class:`Function`
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
        l = test*sqrt(s)*dx
        self.prob = LinearVariationalProblem(a, l, solution)
        self.solver = LinearVariationalSolver(self.prob, solver_parameters=solver_parameters)

    def solve(self):
        """Compute the magnitude"""
        self.solver.solve()
        np.maximum(self.solution.dat.data, self.min_val, self.solution.dat.data)


def compute_bottom_drag(h_b, drag):
    r"""
    Computes bottom drag coefficient (Cd) from the law-of-the wall

    .. math::
        C_D = \left( \frac{\kappa}{\ln (h_b + z_0)/z_0} \right)^2

    :arg h_b: the height above bed where the bottom velocity is evaluated in
        the law-of-the-wall fit
    :type h_b: :class:`Function`
    :arg drag: field where C_D is stored
    :type drag: :class:`Function`
    """
    # FIXME z0 should be a field, i.e. an argument to this function
    von_karman = physical_constants['von_karman']
    z0_friction = physical_constants['z0_friction']
    drag.assign((von_karman / ln((h_b + z0_friction)/z0_friction))**2)
    return drag


def get_horizontal_elem_size_2d(sol2d):
    """
    Computes horizontal element size from the 2D mesh

    :arg sol2d: 2D :class:`Function` where result is stored
    """
    p1_2d = sol2d.function_space()
    mesh = p1_2d.mesh()
    test = TestFunction(p1_2d)
    tri = TrialFunction(p1_2d)
    a = inner(test, tri) * dx
    l = inner(test, sqrt(CellVolume(mesh))) * dx
    solve(a == l, sol2d)


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
    out_function.interpolate(
        Expression('f0+beta*(x[1]-y_0)', f0=f0, beta=beta, y_0=y_offset))

#TODO Edit for 2D only
class SmagorinskyViscosity(object):
    r"""
    Computes Smagorinsky subgrid scale horizontal viscosity

    This formulation is according to Ilicak et al. (2012) and
    Griffies and Hallberg (2000).

    .. math::
        \nu = (C_s \Delta x)^2 |S|

    with the deformation rate

    .. math::
        |S| &= \sqrt{D_T^2 + D_S^2} \\
        D_T &= \frac{\partial u}{\partial x} - \frac{\partial v}{\partial y} \\
        D_S &= \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}

    :math:`\Delta x` is the horizontal element size and :math:`C_s` is the
    Smagorinsky coefficient.

    To match a certain mesh Reynolds number :math:`Re_h` set
    :math:`C_s = 1/\sqrt{Re_h}`.

    Ilicak et al. (2012). Spurious dianeutral mixing and the role of
    momentum closure. Ocean Modelling, 45-46(0):37-58.
    http://dx.doi.org/10.1016/j.ocemod.2011.10.003

    Griffies and Hallberg (2000). Biharmonic friction with a
    Smagorinsky-like viscosity for use in large-scale eddy-permitting
    ocean models. Monthly Weather Review, 128(8):2935-2946.
    http://dx.doi.org/10.1175/1520-0493(2000)128%3C2935:BFWASL%3E2.0.CO;2
    """
    def __init__(self, uv, output, c_s, h_elem_size, max_val, min_val=1e-10,
                 weak_form=True, solver_parameters={}):
        """
        :arg uv_3d: horizontal velocity
        :type uv_3d: 3D vector :class:`Function`
        :arg output: Smagorinsky viscosity field
        :type output: 3D scalar :class:`Function`
        :arg c_s: Smagorinsky coefficient
        :type c_s: float or :class:`Constant`
        :arg h_elem_size: field that defines the horizontal element size
        :type h_elem_size: 3D scalar :class:`Function` or :class:`Constant`
        :arg float max_val: Maximum allowed viscosity. Viscosity will be clipped at
            this value.
        :kwarg float min_val: Minimum allowed viscosity. Viscosity will be clipped at
            this value.
        :kwarg bool weak_form: Compute velocity shear by integrating by parts.
            Necessary for some function spaces (e.g. P0).
        :kwarg dict solver_parameters: PETSc solver options
        """
        solver_parameters.setdefault('ksp_atol', 1e-12)
        solver_parameters.setdefault('ksp_rtol', 1e-16)
        assert max_val.function_space() == output.function_space(), \
            'max_val function must belong to the same space as output'
        self.max_val = max_val
        self.min_val = min_val
        self.output = output
        self.weak_form = weak_form

        if self.weak_form:
            # solve grad(u) weakly
            mesh = output.function_space().mesh()
            fs_grad = FunctionSpace(mesh, 'DP', 1, vfamily='DP', vdegree=1)
            self.grad = {}
            for icomp in [0, 1]:
                for j in [0, 1]:
                    self.grad[(icomp, j)] = Function(fs_grad, name='uv_grad({:},{:})'.format(icomp, j))

            tri_grad = TrialFunction(fs_grad)
            test_grad = TestFunction(fs_grad)

            normal = FacetNormal(mesh)
            a = inner(tri_grad, test_grad)*dx

            self.solver_grad = {}
            for icomp in [0, 1]:
                for j in [0, 1]:
                    a = inner(tri_grad, test_grad)*dx
                    # l = inner(Dx(uv[0], 0), test_grad)*dx
                    l = -inner(Dx(test_grad, j), uv[icomp])*dx
                    l += inner(avg(uv[icomp]), jump(test_grad, normal[j]))*dS_v
                    l += inner(uv[icomp], test_grad*normal[j])*ds_v
                    prob = LinearVariationalProblem(a, l, self.grad[(icomp, j)])
                    self.solver_grad[(icomp, j)] = LinearVariationalSolver(prob, solver_parameters=solver_parameters)

        fs = output.function_space()
        tri = TrialFunction(fs)
        test = TestFunction(fs)

        # rate of strain tensor
        if self.weak_form:
            d_t = self.grad[(0, 0)] - self.grad[(1, 1)]
            d_s = self.grad[(0, 1)] + self.grad[(1, 0)]
        else:
            d_t = Dx(uv[0], 0) - Dx(uv[1], 1)
            d_s = Dx(uv[0], 1) + Dx(uv[1], 0)
        nu = c_s**2*h_elem_size**2 * sqrt(d_t**2 + d_s**2)

        a = test*tri*dx
        l = test*nu*dx
        self.prob = LinearVariationalProblem(a, l, output)
        self.solver = LinearVariationalSolver(self.prob, solver_parameters=solver_parameters)

    def solve(self):
        """Compute viscosity"""
        if self.weak_form:
            for icomp in [0, 1]:
                for j in [0, 1]:
                    self.solver_grad[(icomp, j)].solve()
        self.solver.solve()
        # remove negative values
        ix = self.output.dat.data < self.min_val
        self.output.dat.data[ix] = self.min_val

        # crop too large values
        ix = self.output.dat.data > self.max_val.dat.data
        self.output.dat.data[ix] = self.max_val.dat.data[ix]


def tensor_jump(v, n):
    r"""
    Jump term for vector functions based on the tensor product

    .. math::
        \text{jump}(\mathbf{u}, \mathbf{n}) = (\mathbf{u}^+ \mathbf{n}^+) +
        (\mathbf{u}^- \mathbf{n}^-)

    This is the discrete equivalent of grad(u) as opposed to the
    vectorial UFL jump operator :meth:`ufl.jump` which represents div(u).
    """
    return outer(v('+'), n('+')) + outer(v('-'), n('-'))


def compute_boundary_length(mesh2d):
    """
    Computes the length of the boundary segments in given 2d mesh
    """
    p1 = FunctionSpace(mesh2d, 'CG', 1)
    boundary_markers = mesh2d.exterior_facets.unique_markers
    boundary_len = {}
    for i in boundary_markers:
        ds_restricted = ds(int(i))
        one_func = Function(p1).assign(1.0)
        boundary_len[i] = assemble(one_func * ds_restricted)
    return boundary_len
