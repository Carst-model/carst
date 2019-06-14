import math
import firedrake as fd
from functions import carst_funcs, FunctionContainer
from abc import ABC, abstractmethod
from copy import deepcopy

# Set numerical constants
TINY = 1e-10


class GenericSolver(ABC):
    # Override these attrs to specify what functions and files your solver wants
    _WANTED_FUNCTIONS = ()
    _WANTED_FILES = ()

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
        self.land = land(self.coordinate_space, self.function_space)

        # Initialise funcs and output files
        self._funcs = FunctionContainer(self, type(self)._WANTED_FUNCTIONS)
        self._out_files = {file_name: fd.File("{0}/{1}.pvd".format(self._output_folder, file_name)) for file_name in type(self)._WANTED_FILES}

    def set_condition(self, condition):
        self._funcs[carst_funcs.sed].assign(condition)
        self._funcs[carst_funcs.sed_old].assign(condition)

    @property
    def function_states(self):
        return deepcopy(self._funcs)

    @property
    def time(self):
        return self._current_time

    @abstractmethod
    def advance(self):
        pass

    @abstractmethod
    def output(self):
        pass


class DiffusionSolver(GenericSolver):
    _WANTED_FUNCTIONS = (
        carst_funcs.sed,
        carst_funcs.sed_old,
        carst_funcs.surface,
        carst_funcs.limiter,
        carst_funcs.thickness,
        carst_funcs.depth,
        carst_funcs.diff,
        # carst_funcs.sea_level,
    )
    _WANTED_FILES = ("layer_data", "surfaces", "sea_level", "land")

    def __init__(self, base_mesh, land, sea_level_constant, times, output_folder):
        super().__init__(base_mesh, land, sea_level_constant, times, output_folder)
        self._diffusion_equation = (
            fd.inner(
                (self._funcs[carst_funcs.sed] - self._funcs[carst_funcs.sed_old]) / times[1],
                self.test_function
            ) + self._funcs[carst_funcs.limiter]
            * self._funcs[carst_funcs.diff]
            * fd.inner(
                fd.grad(self._funcs[carst_funcs.sed] + self.land),
                fd.grad(self.test_function)
            )
        ) * fd.dx

        # Add the sea_level function as an interpolation of _sea_level_constant
        self._funcs.add(carst_funcs.sea_level, fd.Function(
            self.function_space,
            val=fd.interpolate(self._sea_level_constant, self.function_space),
            name=str(carst_funcs.sea_level),
        ))

        self._funcs.interpolate(
            carst_funcs.surface,
            carst_funcs.thickness,
            carst_funcs.limiter,
            carst_funcs.depth,
        )
        self._funcs[carst_funcs.diff].project((
            2 / fd.sqrt(2 * math.pi)
            * fd.exp(-0.5 * self._funcs[carst_funcs.depth] ** 2)
        ) + 0.2022)

        # Write out initial data to our blank output files
        self._out_files["land"].write(self.land)
        self.output()

    def advance(self):
        self._current_time += self._time_step

        # Perform the solve
        fd.solve(self._diffusion_equation == 0, self._funcs[carst_funcs.sed])

        # Advance the functions
        self._funcs[carst_funcs.sed_old].assign(self._funcs[carst_funcs.sed])
        self._funcs.interpolate(
            carst_funcs.limiter,
            carst_funcs.surface,
            carst_funcs.depth,
            carst_funcs.diff,
            carst_funcs.thickness,
        )
        self._sea_level_constant.assign(25 * fd.sin(self._current_time / 100000 * 180 / math.pi))
        self._funcs[carst_funcs.sea_level].interpolate(self._sea_level_constant)
        self._funcs[carst_funcs.slf].interpolate(self._funcs[carst_funcs.sea_level])

        # If the current cycle is one where we should output, write to our files
        if self._current_time % self._output_time == 0:
            self.output()

    def output(self):
        self._out_files["layer_data"].write(
            self._funcs[carst_funcs.diff],
            self._funcs[carst_funcs.thickness],
            self._funcs[carst_funcs.depth],
            time=self._current_time,
        )
        self._out_files["surfaces"].write(
            self._funcs[carst_funcs.surface],
            self._funcs[carst_funcs.sed],
            time=self._current_time,
        )
        self._out_files["sea_level"].write(
            self._funcs[carst_funcs.slf],
            time=self._current_time,
        )


# class CarbonateSolver(GenericSolver):
    # _wanted_functions = (
        # carst_funcs.sed,
        # carst_funcs.sed_old,
        # carst_funcs.surface,
        # carst_funcs.limiter,
        # carst_funcs.thickness,
        # carst_funcs.depth,
        # carst_funcs.diff,
        # carst_funcs.light_attenuation,
    # )
    # _wanted_files = ("layer_data", "sea_level", "land")

    # def __init__(self, base_mesh, land, output_folder):
        # super().__init__(base_mesh, land, output_folder)

    # def diffuse_carbonate_test(self, initial_condition, start_time, end_time, time_step, output_time, carb_pr):
        # # Ensure we have a land attr set
        # if not hasattr(self, "land"):
            # raise AttributeError("No land passed to module")

        # # Initialise our functions and out files
        # funcs = FunctionContainer(self, DiffusionSolver._wanted_functions["carbonates"])
        # out_files = self._get_files("carbonates")
        # sl_time_constant = fd.Constant(25 * fd.sin(start_time / 1000000 / 180 * math.pi))
        # funcs[carst_funcs.sea_level] = fd.Function(
            # self.test_function,
            # val=fd.interpolate(sl_time_constant, self.test_function),
            # name="sea_level",
        # )
        # funcs.interpolate(carst_funcs.sea_level)

        # out_files["land"].write(self.land)

        # # Assign our two sediment buildup functions to the initial condition
        # funcs[carst_funcs.sed].assign(initial_condition)
        # funcs[carst_funcs.sed_old].assign(initial_condition)

        # # Initial advancement
        # funcs.interpolate((
            # carst_funcs.surface,
            # carst_funcs.thickness,
            # carst_funcs.limiter,
            # carst_funcs.depth,
        # ))
        # funcs[carst_funcs.diff].project((
            # 2 / fd.sqrt(2 * math.pi)
            # * fd.exp(-0.5 * funcs[carst_funcs.depth] ** 2)
        # ) + 0.2022)

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

        # # Write out initial data to our blank output files
        # out_files["surfaces"].write(funcs[carst_funcs.surface], funcs[carst_funcs.sed], time=0)
        # out_files["layer_data"].write(funcs[carst_funcs.diff], funcs[carst_funcs.thickness], funcs[carst_funcs.depth], time=0)
        # out_files["sea_level"].write(funcs[carst_funcs.slf], time=0)

        # # Begin looping until our end time has been reached
        # current_time = start_time
        # while current_time <= end_time:
            # current_time += time_step

            # # Perform the solve
            # funcs.interpolate(carst_funcs.light_attenuation)
            # fd.solve(F == 0, funcs[carst_funcs.sed])
            # funcs[carst_funcs.sed] += carb_pr * funcs[carst_funcs.light_attenuation]

            # # Iterate
            # funcs[carst_funcs.sed_old].assign(funcs[carst_funcs.sed])
            # # Advance the functions
            # funcs.interpolate((
                # carst_funcs.limiter,
                # carst_funcs.surface,
                # carst_funcs.depth,
                # carst_funcs.diff,
                # carst_funcs.thickness,
            # ))
            # sl_time_constant.assign(25 * fd.sin(current_time / 100000 * 180 / math.pi))
            # funcs[carst_funcs.sea_level].interpolate(sl_time_constant)
            # funcs[carst_funcs.slf].interpolate(funcs[carst_funcs.sea_level])

            # # If the current cycle is one where we should output, write to our files
            # if current_time % output_time == 0:
                # out_files["layer_data"].write(
                    # funcs[carst_funcs.diff],
                    # funcs[carst_funcs.thickness],
                    # funcs[carst_funcs.depth],
                    # time=current_time
                # )
                # out_files["surfaces"].write(
                    # funcs[carst_funcs.surface],
                    # funcs[carst_funcs.sed],
                    # time=current_time
                # )
                # out_files["sea_level"].write(
                    # funcs[carst_funcs.slf],
                    # time=current_time
                # )
