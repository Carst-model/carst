import math
import firedrake as fd
from functions import carst_funcs, FunctionContainer
from abc import ABC, abstractmethod
from copy import deepcopy

# Set numerical constants
TINY = 1e-10


class GenericSolver(ABC):
    # Override these attrs to specify what your solver needs
    # Functions the solver needs to have
    WANTED_FUNCTIONS = ()
    # File the solver wants to write to (along with a function to return a list of things to write to it
    _WANTED_FILES = {
        "land": None,
    }
    # The order the functions should be interpolated when the class is instantiated
    _INIT_INTERPOLATION_ORDER = ()
    # The order the functions should be interpolated in advance()
    _INTERPOLATION_ORDER = ()

    def __init__(self, base_mesh, land, sea_level_constant, times, output_folder):
        # Need type checking here
        # Need to check whether the output folder actually exists with os
        # Store the passed values
        self._sea_level_constant = sea_level_constant
        self._current_time, self._time_step, self._output_time = times
        self._output_folder = output_folder
        self.mesh = base_mesh

        # Generate our workspace from the mesh
        self.coordinate_space = fd.SpatialCoordinate(self.mesh)
        self.function_space = fd.FunctionSpace(self.mesh, "CG", 1)
        self.test_function = fd.TestFunction(self.function_space)

        # Execute the land function
        self._land = land(self.coordinate_space, self.function_space)

        # Initialise funcs and perform the first iteration
        self._funcs = FunctionContainer(self)
        self._funcs.interpolate(*type(self)._INIT_INTERPOLATION_ORDER)

        # Initialise _out_files and write land to them
        self._out_files = {file_name: fd.File("{0}/{1}.pvd".format(self._output_folder, file_name)) for file_name in type(self)._WANTED_FILES}
        self._out_files["land"].write(self.land)

    def set_condition(self, condition):
        self._funcs[carst_funcs.sed].assign(condition)
        self._funcs[carst_funcs.sed_old].assign(condition)

    @property
    def function_states(self):
        return deepcopy(self._funcs)

    @property
    def time(self):
        return self._current_time

    @property
    def time_step(self):
        return self._time_step

    @property
    def sea_level_constant(self):
        # return deepcopy(self._sea_level_constant)
        return self._sea_level_constant

    @property
    def land(self):
        return self.land
        # return deepcopy(self.land)

    @abstractmethod
    def advance(self):
        pass

    def output(self):
        for file_name in type(self)._WANTED_FILES.keys():
            if file_name != "land":
                self._out_files[file_name].write(
                    *type(self)._WANTED_FILES[file_name](self._funcs),
                    time=self._current_time,
                )


class DiffusionSolver(GenericSolver):
    WANTED_FUNCTIONS = super().WANTED_FUNCTIONS + (
        carst_funcs.sed,
        carst_funcs.sed_old,
        carst_funcs.surface,
        carst_funcs.limiter,
        carst_funcs.thickness,
        carst_funcs.depth,
        carst_funcs.diff,
        carst_funcs.sea_level,
    )
    _INIT_INTERPOLATION_ORDER = super()._INIT_INTERPOLATION_ORDER + (
        carst_funcs.surface,
        carst_funcs.thickness,
        carst_funcs.limiter,
        carst_funcs.depth,
        carst_funcs.sea_level,
    )
    _INTERPOLATION_ORDER = super()._INTERPOLATION_ORDER + (
        carst_funcs.limiter,
        carst_funcs.surface,
        carst_funcs.depth,
        carst_funcs.diff,
        carst_funcs.thickness,
    )

    _WANTED_FILES = super()._WANTED_FILES + {
        "layer_data": lambda funcs: (funcs[carst_funcs.diff], funcs[carst_funcs.thickness], funcs[carst_funcs.depth]),
        "surfaces": lambda funcs: (funcs[carst_funcs.surface], funcs[carst_funcs.sed]),
        "sea_level": lambda funcs: (funcs[carst_funcs.sea_level]),
    }

    _DIFFUSION_EQUATION_GENERIC = lambda solver, funcs: (
        fd.inner(
            (funcs[carst_funcs.sed] - funcs[carst_funcs.sed_old]) / solver.time_step,
            solver.test_function
        ) + funcs[carst_funcs.limiter]
        * funcs[carst_funcs.diff]
        * fd.inner(
            fd.grad(funcs[carst_funcs.sed] + solver.land),
            fd.grad(solver.test_function)
        )
    ) * fd.dx

    def __init__(self, base_mesh, land, sea_level_constant, times, output_folder):
        super().__init__(base_mesh, land, sea_level_constant, times, output_folder)
        self._diffusion_equation = type(self)._DIFFUSION_EQUATION_GENERIC(self, self._funcs)

        # Add the sea_level function as an interpolation of _sea_level_constant
        self._funcs[carst_funcs.sea_level].interpolate(self._sea_level_constant)

        self._funcs[carst_funcs.diff].project((
            2 / fd.sqrt(2 * math.pi)
            * fd.exp(-0.5 * self._funcs[carst_funcs.depth] ** 2)
        ) + 0.2022)

        # Write out initial data to our blank output files
        self.output()

    def advance(self):
        self._current_time += self._time_step

        # Perform the solve
        fd.solve(self._diffusion_equation == 0, self._funcs[carst_funcs.sed])

        # Advance the functions
        self._funcs[carst_funcs.sed_old].assign(self._funcs[carst_funcs.sed])
        self._funcs.interpolate(*DiffusionSolver._INTERPOLATION_ORDER)
        self._sea_level_constant.assign(25 * fd.sin(self._current_time / 100000 * 180 / math.pi))
        self._funcs[carst_funcs.sea_level].interpolate(self._sea_level_constant)

        # If the current cycle is one where we should output, write to our files
        if self._current_time % self._output_time == 0:
            self.output()


class CarbonateSolver(GenericSolver):
    WANTED_FUNCTIONS = super().WANTED_FUNCTIONS + (
        carst_funcs.sed,
        carst_funcs.sed_old,
        carst_funcs.surface,
        carst_funcs.limiter,
        carst_funcs.thickness,
        carst_funcs.depth,
        carst_funcs.diff,
        carst_funcs.light_attenuation,
        carst_funcs.sea_level,
    )
    _INTERPOLATION_ORDER = super()._INTERPOLATION_ORDER + (
        carst_funcs.limiter,
        carst_funcs.surface,
        carst_funcs.depth,
        carst_funcs.diff,
        carst_funcs.thickness,
    )
    _INIT_INTERPOLATION_ORDER = super()._INIT_INTERPOLATION_ORDER + (
        carst_funcs.surface,
        carst_funcs.thickness,
        carst_funcs.limiter,
        carst_funcs.depth,
    )
    _WANTED_FILES = super()._WANTED_FILES + ("layer_data", "surfaces", "sea_level")

    def __init__(self, base_mesh, land, sea_level_constant, times, output_folder, carb_pr):
        super().__init__(base_mesh, land, sea_level_constant, times, output_folder)

        self._carb_pr = carb_pr

        # Initial advancement
        self._funcs[carst_funcs.sea_level].interpolate(self._sea_level_constant)
        self._funcs[carst_funcs.diff].project((
            2 / fd.sqrt(2 * math.pi)
            * fd.exp(-0.5 * self._funcs[carst_funcs.depth] ** 2)
        ) + 0.2022)
        self.output()

    def advance(self):
        # F = (
            # fd.inner(
                # (funcs[carst_funcs.sed] - funcs[carst_funcs.sed_old]) / time_step,
                # self.test_function
            # ) + funcs[carst_funcs.limiter]
            # * funcs[carst_funcs.diff]
            # * fd.inner(
                # fd.grad(funcs[carst_funcs.sed] + self.land),
                # fd.grad(self.test_function)
            # )
        # ) * fd.dx

        # Begin looping until our end time has been reached
        self._current_time += self._time_step

        # Perform the solve
        self._funcs.interpolate(carst_funcs.light_attenuation)
        # fd.solve(F == 0, funcs[carst_funcs.sed])
        self._funcs[carst_funcs.sed] += self._carb_pr * self._funcs[carst_funcs.light_attenuation]

        # Iterate
        self._funcs[carst_funcs.sed_old].assign(self._funcs[carst_funcs.sed])
        # Advance the functions
        self._funcs.interpolate(*type(self)._INTERPOLATION_ORDER)
        self._sea_level_constant.assign(25 * fd.sin(self._current_time / 100000 * 180 / math.pi))
        self._funcs[carst_funcs.sea_level].interpolate(self._sea_level_constant)

        # If the current cycle is one where we should output, write to our files
        if self._current_time % self._output_time == 0:
            self.output()
