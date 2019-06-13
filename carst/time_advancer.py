import math
import firedrake as fd
import os
from functions import carst_funcs, FunctionContainer

# Set numerical constants
TINY = 1e-10


class DiffusionSolver:
    # Declare the names of the functions we want for each solver
    _wanted_functions = {
        "real_scale": (
            carst_funcs.sed,
            carst_funcs.sed_old,
            carst_funcs.surface,
            carst_funcs.limiter,
            carst_funcs.thickness,
            carst_funcs.depth,
            carst_funcs.diff,
        ),
        "carbonates": (
            carst_funcs.sed,
            carst_funcs.sed_old,
            carst_funcs.surface,
            carst_funcs.limiter,
            carst_funcs.thickness,
            carst_funcs.depth,
            carst_funcs.diff,
            carst_funcs.light_attenuation,
        ),
    }

    # Declare the names of the files the solvers each want
    _wanted_files = {
        "real_scale": ("surfaces", "layer_data", "sea_level", "land"),
        "carbonates": ("surfaces", "layer_data", "sea_level", "land"),
    }

    def __init__(self, base_mesh, land, output_folder):
        # Need to check whether the output folder actually exists with os
        self.output_folder = output_folder

        # Generate our workspace from the mesh
        self.mesh = base_mesh
        self.coordinate_space = fd.SpatialCoordinate(self.mesh)
        self.function_space = fd.FunctionSpace(self.mesh, "CG", 1)
        self.test_function = fd.TestFunction(self.function_space)

        self.land = land(self.coordinate_space, self.function_space)

    def _get_files(self, func_group):
        return {file_name: fd.File("{0}/{1}.pvd".format(self.output_folder, file_name)) for file_name in DiffusionSolver._wanted_files[func_group]}

    def diffuse_real_scale_test(self, initial_condition, start_time, end_time, time_step, output_time):
        if not hasattr(self, "land"):
            raise AttributeError("No land passed to module")

        sl_time_constant = fd.Constant(25 * fd.sin(start_time / 100000 / 180 * math.pi))

        # Initialise our functions and out files
        out_files = self._get_files("real_scale")
        funcs = FunctionContainer(self, DiffusionSolver._wanted_functions["real_scale"])

        # Currently not easy to slot this workaround into FunctionContainer so it can stay here for now
        # if we don't do this, the sea level ouput has a random function name.
        # Not helpful...
        funcs.functions[carst_funcs.slf] = fd.Function(self.function_space, name="sea_level")
        funcs.functions[carst_funcs.sea_level] = fd.Function(
            self.function_space,
            val=fd.interpolate(sl_time_constant, self.function_space),
            name="sea_level",
        )
        funcs.functions[carst_funcs.slf].interpolate(funcs.functions[carst_funcs.sea_level])

        out_files["land"].write(self.land)

        # Assign our two sediment buildup functions to the initial condition
        funcs.functions[carst_funcs.sed].assign(initial_condition)
        funcs.functions[carst_funcs.sed_old].assign(initial_condition)

        funcs.interpolate((
            carst_funcs.surface,
            carst_funcs.thickness,
            carst_funcs.limiter,
            carst_funcs.depth,
        ))
        funcs.functions[carst_funcs.diff].project((
            2 / fd.sqrt(2 * math.pi)
            * fd.exp(-0.5 * funcs.functions[carst_funcs.depth] ** 2)
        ) + 0.2022)

        # Specify our function to solve (ie. the diffusion equation)
        F = (
            fd.inner(
                (funcs.functions[carst_funcs.sed] - funcs.functions[carst_funcs.sed_old]) / time_step,
                self.test_function
            ) + funcs.functions[carst_funcs.limiter]
            * funcs.functions[carst_funcs.diff]
            * fd.inner(
                fd.grad(funcs.functions[carst_funcs.sed] + self.land),
                fd.grad(self.test_function)
            )
        ) * fd.dx

        # Write out initial data to our blank output files
        out_files["surfaces"].write(funcs.functions[carst_funcs.surface], funcs.functions[carst_funcs.sed], time=0)
        out_files["layer_data"].write(funcs.functions[carst_funcs.diff], funcs.functions[carst_funcs.thickness], funcs.functions[carst_funcs.depth], time=0)
        out_files["sea_level"].write(funcs.functions[carst_funcs.slf], time=0)

        # Begin looping until our end time has been reached
        current_time = start_time
        while current_time <= end_time:
            current_time += time_step

            # Perform the solve
            fd.solve(F == 0, funcs.functions[carst_funcs.sed])

            # Iterate
            funcs.functions[carst_funcs.sed_old].assign(funcs.functions[carst_funcs.sed])
            # Advance the functions
            funcs.interpolate((
                carst_funcs.limiter,
                carst_funcs.surface,
                carst_funcs.depth,
                carst_funcs.diff,
                carst_funcs.thickness,
            ))
            sl_time_constant.assign(25 * fd.sin(current_time / 100000 * 180 / math.pi))
            funcs.functions[carst_funcs.sea_level].interpolate(sl_time_constant)
            funcs.functions[carst_funcs.slf].interpolate(funcs.functions[carst_funcs.sea_level])

            # If the current cycle is one where we should output, write to our files
            if current_time % output_time == 0:
                out_files["layer_data"].write(
                    funcs.functions[carst_funcs.diff],
                    funcs.functions[carst_funcs.thickness],
                    funcs.functions[carst_funcs.depth],
                    time=current_time
                )
                out_files["surfaces"].write(
                    funcs.functions[carst_funcs.surface],
                    funcs.functions[carst_funcs.sed],
                    time=current_time
                )
                out_files["sea_level"].write(
                    funcs.functions[carst_funcs.slf],
                    time=current_time
                )

    def diffuse_carbonate_test(self, initial_condition, start_time, end_time, time_step, output_time, carb_pr):
        # Ensure we have a land attr set
        if not hasattr(self, "land"):
            raise AttributeError("No land passed to module")

        # Initialise our functions and out files
        funcs = FunctionContainer(self, DiffusionSolver._wanted_functions["carbonates"])
        out_files = self._get_files("carbonates")
        sl_time_constant = fd.Constant(25 * fd.sin(start_time / 1000000 / 180 * math.pi))
        funcs.functions[carst_funcs.sea_level] = fd.Function(
            self.test_function,
            val=fd.interpolate(sl_time_constant, self.test_function),
            name="sea_level",
        )
        funcs.interpolate(carst_funcs.sea_level)

        out_files["land"].write(self.land)

        # Assign our two sediment buildup functions to the initial condition
        funcs.functions[carst_funcs.sed].assign(initial_condition)
        funcs.functions[carst_funcs.sed_old].assign(initial_condition)

        # Initial advancement
        funcs.interpolate((
            carst_funcs.surface,
            carst_funcs.thickness,
            carst_funcs.limiter,
            carst_funcs.depth,
        ))
        funcs.functions[carst_funcs.diff].project((
            2 / fd.sqrt(2 * math.pi)
            * fd.exp(-0.5 * funcs.functions[carst_funcs.depth] ** 2)
        ) + 0.2022)

        F = (
            fd.inner(
                (funcs.functions[carst_funcs.sed] - funcs.functions[carst_funcs.sed_old]) / time_step,
                self.test_function
            ) + funcs.functions[carst_funcs.limiter]
            * funcs.functions[carst_funcs.diff]
            * fd.inner(
                fd.grad(funcs.functions[carst_funcs.sed] + self.land),
                fd.grad(self.test_function)
            )
        ) * fd.dx

        # Write out initial data to our blank output files
        out_files["surfaces"].write(funcs.functions[carst_funcs.surface], funcs.functions[carst_funcs.sed], time=0)
        out_files["layer_data"].write(funcs.functions[carst_funcs.diff], funcs.functions[carst_funcs.thickness], funcs.functions[carst_funcs.depth], time=0)
        out_files["sea_level"].write(funcs.functions[carst_funcs.slf], time=0)

        # Begin looping until our end time has been reached
        current_time = start_time
        while current_time <= end_time:
            current_time += time_step

            # Perform the solve
            funcs.interpolate(carst_funcs.light_attenuation)
            fd.solve(F == 0, funcs.functions[carst_funcs.sed])
            funcs.functions[carst_funcs.sed] += carb_pr * funcs.functions[carst_funcs.light_attenuation]

            # Iterate
            funcs.functions[carst_funcs.sed_old].assign(funcs.functions[carst_funcs.sed])
            # Advance the functions
            funcs.interpolate((
                carst_funcs.limiter,
                carst_funcs.surface,
                carst_funcs.depth,
                carst_funcs.diff,
                carst_funcs.thickness,
            ))
            sl_time_constant.assign(25 * fd.sin(current_time / 100000 * 180 / math.pi))
            funcs.functions[carst_funcs.sea_level].interpolate(sl_time_constant)
            funcs.functions[carst_funcs.slf].interpolate(funcs.functions[carst_funcs.sea_level])

            # If the current cycle is one where we should output, write to our files
            if current_time % output_time == 0:
                out_files["layer_data"].write(
                    funcs.functions[carst_funcs.diff],
                    funcs.functions[carst_funcs.thickness],
                    funcs.functions[carst_funcs.depth],
                    time=current_time
                )
                out_files["surfaces"].write(
                    funcs.functions[carst_funcs.surface],
                    funcs.functions[carst_funcs.sed],
                    time=current_time
                )
                out_files["sea_level"].write(
                    funcs.functions[carst_funcs.slf],
                    time=current_time
                )
