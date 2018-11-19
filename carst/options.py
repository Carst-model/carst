"""
Thetis options for the 2D and 3D model

All options are type-checked and they are stored in traitlets Configurable
objects.
"""
from .configuration import *
from .firedrake import Constant


class TimeStepperOptions(FrozenHasTraits):
    """Base class for all time stepper options"""
    name = 'Time stepper'


class ExplicitTimestepperOptions(TimeStepperOptions):
    """Options for explicit time integrator"""
    use_automatic_timestep = Bool(True, help='Set time step automatically based on local CFL conditions.').tag(config=True)


class SemiImplicitTimestepperOptions2d(TimeStepperOptions):
    """Options for 2d explicit time integrator"""
    solver_parameters = PETScSolverParameters({
        'ksp_type': 'gmres',
        'pc_type': 'fieldsplit',
        'pc_fieldsplit_type': 'multiplicative',
    }).tag(config=True)
    solver_parameters_tracer = PETScSolverParameters({
        'ksp_type': 'gmres',
        'pc_type': 'sor',
    }).tag(config=True)


class SteadyStateTimestepperOptions2d(TimeStepperOptions):
    """Options for 2d steady state solver"""
    solver_parameters = PETScSolverParameters({
        'ksp_type': 'preonly',
        'pc_type': 'lu',
        'mat_type': 'aij'
    }).tag(config=True)


class CrankNicolsonTimestepperOptions2d(SemiImplicitTimestepperOptions2d):
    """Options for 2d Crank-Nicolson time integrator"""
    implicitness_theta = BoundedFloat(
        default_value=0.5, bounds=[0.5, 1.0],
        help='implicitness parameter theta. Value 0.5 implies Crank-Nicolson scheme, 1.0 implies fully implicit formulation.').tag(config=True)
    use_semi_implicit_linearization = Bool(
        False, help="Use linearized semi-implicit time integration").tag(config=True)


class PressureProjectionTimestepperOptions2d(TimeStepperOptions):
    """Options for 2d pressure-projection time integrator"""
    solver_parameters_pressure = PETScSolverParameters({
        'ksp_type': 'preonly',  # we solve the full schur complement exactly, so no need for outer krylov
        'mat_type': 'matfree',
        'pc_type': 'fieldsplit',
        'pc_fieldsplit_type': 'schur',
        'pc_fieldsplit_schur_fact_type': 'full',
        # velocity mass block:
        'fieldsplit_U_2d': {
            'ksp_type': 'gmres',
            'pc_type': 'python',
            'pc_python_type': 'firedrake.AssembledPC',
            'assembled_ksp_type': 'preonly',
            'assembled_pc_type': 'bjacobi',
            'assembled_sub_pc_type': 'ilu',
        },
        # schur system: explicitly assemble the schur system
        # this only works with pressureprojectionicard if the velocity block is just the mass matrix
        # and if the velocity is DG so that this mass matrix can be inverted explicitly
        'fieldsplit_H_2d': {
            'ksp_type': 'preonly',
            'pc_type': 'python',
            'pc_python_type': 'thetis.AssembledSchurPC',
            'schur_ksp_type': 'gmres',
            'schur_ksp_max_it': 100,
            'schur_ksp_converged_reason': True,
            'schur_pc_type': 'gamg',
        },
    }).tag(config=True)
    solver_parameters_momentum = PETScSolverParameters({
        'ksp_type': 'gmres',
        'ksp_converged_reason': True,
        'pc_type': 'bjacobi',
        'sub_ksp_type': 'preonly',
        'sub_pc_type': 'sor',
    }).tag(config=True)
    implicitness_theta = BoundedFloat(
        default_value=0.5, bounds=[0.5, 1.0],
        help='implicitness parameter theta. Value 0.5 implies Crank-Nicolson scheme, 1.0 implies fully implicit formulation.').tag(config=True)
    use_semi_implicit_linearization = Bool(
        True, help="Use linearized semi-implicit time integration").tag(config=True)
    picard_iterations = PositiveInteger(
        default_value=2,
        help='number of Picard iterations to converge the nonlinearity in the equations.')


class ExplicitTimestepperOptions2d(ExplicitTimestepperOptions):
    """Options for 2d explicit time integrator"""
    solver_parameters = PETScSolverParameters({
        'snes_type': 'ksponly',
        'ksp_type': 'cg',
        'pc_type': 'bjacobi',
        'sub_ksp_type': 'preonly',
        'sub_pc_type': 'ilu',
        'mat_type': 'aij',
    }).tag(config=True)



class CommonModelOptions(FrozenConfigurable):
    """Options that are common"""
    name = 'Model options'
    polynomial_degree = NonNegativeInteger(1, help='Polynomial degree of elements').tag(config=True)
    element_family = Enum(
        ['dg-dg', 'rt-dg', 'dg-cg'],
        default_value='dg-dg',
        help="""Finite element family

        2D solver supports 'dg-dg', 'rt-dg', or 'dg-cg' velocity-pressure pairs.""").tag(config=True)

    use_nonlinear_equations = Bool(True, help='Use nonlinear shallow water equations').tag(config=True)
    use_grad_div_viscosity_term = Bool(
        False,
        help=r"""Include :math:`\nabla (\nu_h \nabla \cdot \bar{\textbf{u}})` term in the depth-averaged viscosity

        See :class:`.shallowwater_eq.HorizontalViscosityTerm` for details.""").tag(config=True)
    use_grad_depth_viscosity_term = Bool(
        True,
        help=r"""Include :math:`\nabla H` term in the depth-averaged viscosity

        See :class:`.shallowwater_eq.HorizontalViscosityTerm` for details.""").tag(config=True)

    use_lax_friedrichs_velocity = Bool(
        True, help="use Lax Friedrichs stabilisation in horizontal momentum advection.").tag(config=True)
    lax_friedrichs_velocity_scaling_factor = FiredrakeConstantTraitlet(
        Constant(1.0), help="Scaling factor for Lax Friedrichs stabilisation term in horiozonal momentum advection.").tag(config=True)
    use_lax_friedrichs_tracer = Bool(
        True, help="Use Lax Friedrichs stabilisation in tracer advection.").tag(config=True)
    lax_friedrichs_tracer_scaling_factor = FiredrakeConstantTraitlet(
        Constant(1.0), help="Scaling factor for tracer Lax Friedrichs stability term.").tag(config=True)
    use_limiter_for_tracers = Bool(
        True, help="Apply P1DG limiter for tracer fields").tag(config=True)

    check_volume_conservation_2d = Bool(
        False, help="""
        Compute volume of the 2D mode at every export

        2D volume is defined as the integral of the water elevation field.
        Prints deviation from the initial volume to stdout.
        """).tag(config=True)
    log_output = Bool(
        True, help="Redirect all output to log file in output directory").tag(config=True)
    timestep = PositiveFloat(
        10.0, help="Time step").tag(config=True)
    cfl_2d = PositiveFloat(
        1.0, help="Factor to scale the 2d time step OBSOLETE").tag(config=True)  # TODO OBSOLETE
    cfl_3d = PositiveFloat(
        1.0, help="Factor to scale the 2d time step OBSOLETE").tag(config=True)  # TODO OBSOLETE
    simulation_export_time = PositiveFloat(
        100.0, help="""
        Export interval in seconds

        All fields in fields_to_export list will be stored to disk and
        diagnostics will be computed
        """).tag(config=True)
    simulation_end_time = PositiveFloat(
        1000.0, help="Simulation duration in seconds").tag(config=True)
    horizontal_velocity_scale = FiredrakeConstantTraitlet(
        Constant(0.1), help="""
        Maximum horizontal velocity magnitude

        Used to compute max stable advection time step.
        """).tag(config=True)
    horizontal_viscosity_scale = FiredrakeConstantTraitlet(
        Constant(1.0), help="""
        Maximum horizontal viscosity

        Used to compute max stable diffusion time step.
        """).tag(config=True)
    output_directory = Unicode(
        'outputs', help="Directory where model output files are stored").tag(config=True)
    no_exports = Bool(
        False, help="""
        Do not store any outputs to disk

        Disables VTK and HDF5 field outputs. and HDF5 diagnostic outputs.
        Used in CI test suite.
        """).tag(config=True)
    export_diagnostics = Bool(
        True, help="Store diagnostic variables to disk in HDF5 format").tag(config=True)
    fields_to_export = List(
        trait=Unicode,
        default_value=['elev_2d', 'uv_2d', 'uv_3d', 'w_3d'],
        help="Fields to export in VTK format").tag(config=True)
    fields_to_export_hdf5 = List(
        trait=Unicode,
        default_value=[],
        help="Fields to export in HDF5 format").tag(config=True)
    verbose = Integer(0, help="Verbosity level").tag(config=True)
    linear_drag_coefficient = FiredrakeScalarExpression(
        None, allow_none=True, help=r"""
        2D linear drag parameter :math:`L`

        Bottom stress is :math:`\tau_b/\rho_0 = -L \mathbf{u} H`
        """).tag(config=True)
    quadratic_drag_coefficient = FiredrakeScalarExpression(
        None, allow_none=True, help=r"""
        Dimensionless 2D quadratic drag parameter :math:`C_D`

        Bottom stress is :math:`\tau_b/\rho_0 = -C_D |\mathbf{u}|\mathbf{u}`
        """).tag(config=True)
    manning_drag_coefficient = FiredrakeScalarExpression(
        None, allow_none=True, help=r"""
        Manning-Strickler 2D quadratic drag parameter :math:`\mu`

        Bottom stress is :math:`\tau_b/\rho_0 = -g \mu^2 |\mathbf{u}|\mathbf{u}/H^{1/3}`
        """).tag(config=True)
    horizontal_viscosity = FiredrakeScalarExpression(
        None, allow_none=True, help="Horizontal viscosity").tag(config=True)
    coriolis_frequency = FiredrakeScalarExpression(
        None, allow_none=True, help="2D Coriolis parameter").tag(config=True)
    wind_stress = FiredrakeVectorExpression(
        None, allow_none=True, help="Stress at free surface (2D vector function)").tag(config=True)
    atmospheric_pressure = FiredrakeScalarExpression(
        None, allow_none=True, help="Atmospheric pressure at free surface, in pascals").tag(config=True)
    momentum_source_2d = FiredrakeVectorExpression(
        None, allow_none=True, help="Source term for 2D momentum equation").tag(config=True)
    volume_source_2d = FiredrakeScalarExpression(
        None, allow_none=True, help="Source term for 2D continuity equation").tag(config=True)
    tracer_source_2d = FiredrakeScalarExpression(
        None, allow_none=True, help="Source term for 2D tracer equation").tag(config=True)
    horizontal_diffusivity = FiredrakeCoefficient(
        None, allow_none=True, help="Horizontal diffusivity for tracers").tag(config=True)


# NOTE all parameters are now case sensitive
# TODO rename time stepper types? Allow capitals and spaces?
@attach_paired_options("timestepper_type",
                       PairedEnum([('SSPRK33', ExplicitTimestepperOptions2d),
                                   ('ForwardEuler', ExplicitTimestepperOptions2d),
                                   ('BackwardEuler', SemiImplicitTimestepperOptions2d),
                                   ('CrankNicolson', CrankNicolsonTimestepperOptions2d),
                                   ('DIRK22', SemiImplicitTimestepperOptions2d),
                                   ('DIRK33', SemiImplicitTimestepperOptions2d),
                                   ('SteadyState', SteadyStateTimestepperOptions2d),
                                   ('PressureProjectionPicard', PressureProjectionTimestepperOptions2d),
                                   ('SSPIMEX', SemiImplicitTimestepperOptions2d),
                                   ],
                                  "timestepper_options",
                                  default_value='CrankNicolson',
                                  help='Name of the time integrator').tag(config=True),
                       Instance(TimeStepperOptions, args=()).tag(config=True))
class ModelOptions2d(CommonModelOptions):
    """Options for 2D depth-averaged shallow water model"""
    name = 'Depth-averaged 2D model'
    solve_tracer = Bool(False, help='Solve tracer transport').tag(config=True)
    use_wetting_and_drying = Bool(
        False, help=r"""bool: Turn on wetting and drying

        Uses the wetting and drying scheme from Karna et al (2011).
        If ``True``, one should also set :attr:`wetting_and_drying_alpha` to control the bathymetry displacement.
        """).tag(config=True)
    wetting_and_drying_alpha = FiredrakeConstantTraitlet(
        Constant(0.5), help=r"""
        Coefficient: Wetting and drying parameter :math:`\alpha`.

        Used in bathymetry displacement function that ensures positive water depths. Unit is meters.
        """).tag(config=True)
    tidal_turbine_farms = Dict(trait=TidalTurbineFarmOptions,
                               default_value={}, help='Dictionary mapping subdomain ids to the options of the corresponding farm')

    check_tracer_conservation = Bool(
        False, help="""
        Compute total tracer mass at every export

        Prints deviation from the initial mass to stdout.
        """).tag(config=True)

    check_tracer_overshoot = Bool(
        False, help="""
        Compute tracer overshoots at every export

        Prints overshoot values that exceed the initial range to stdout.
        """).tag(config=True)
    tracer_only = Bool(
        False, help="""Hold shallow water variables in initial state

        Advects tracer in the associated (constant) velocity field.
        """).tag(config=True)


@attach_paired_options("timestepper_type",
                       PairedEnum([('LeapFrog', ExplicitTimestepperOptions3d),
                                   ('SSPRK22', ExplicitTimestepperOptions3d),
                                   ],
                                  "timestepper_options",
                                  default_value='SSPRK22',
                                  help='Name of the time integrator').tag(config=True),
                       Instance(TimeStepperOptions, args=()).tag(config=True))
@attach_paired_options("turbulence_model_type",
                       PairedEnum([('gls', GLSModelOptions),
                                   ('pacanowski', PacanowskiPhilanderModelOptions)
                                   ],
                                  "turbulence_model_options",
                                  default_value='gls',
                                  help='Type of vertical turbulence model').tag(config=True),
                       Instance(TurbulenceModelOptions, args=()).tag(config=True))
@attach_paired_options("equation_of_state_type",
                       PairedEnum([('full', EquationOfStateOptions),
                                   ('linear', LinearEquationOfStateOptions)],
                                  "equation_of_state_options",
                                  default_value='full',
                                  help='Type of equation of state').tag(config=True),
                       Instance(EquationOfStateOptions, args=()).tag(config=True))
