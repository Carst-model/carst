import math
import firedrake as fd
import os

# Set numerical constants
TINY = 1e-10


class Diffuse_Solver:
    def __init__(self, base_mesh, output_folder):
        # Check whether this exists please
        self.output_folder = output_folder

        # Generate all the needed stuff from the mesh
        self.mesh = base_mesh
        self.coordinate_space = fd.SpatialCoordinate(self.mesh)
        self.function_space = fd.FunctionSpace(self.mesh, "CG", 1)
        self.test_function = fd.TestFunction(self.function_space)

        # Declare the names of the functions we want for each solver
        self._wanted_functions = {
            "real_scale": (
                "sed",
                "sed_old",
                "surface",
                "limiter",
                "thickness",
                "depth",
                "diff",
            ),
            "carbonates": (
                "sed",
                "sed_old",
                "surface",
                "limiter",
                "thickness",
                "depth",
                "diff",
                "light_attenuation",
            ),
        }

        # Declare the names of the files the solvers each want
        self._wanted_files = {
            "real_scale": ("surfaces", "layer_data", "sea_level", "land"),
            "carbonates": ("surfaces", "layer_data", "sea_level", "land"),
        }

    def add_land(self, land):
        self.land = land

    def _initialise_features(self, func_group):
        return (
            {func_name: fd.Function(self.function_space, name=func_name) for func_name in self._wanted_functions[func_group]},
            {file_name: fd.File("{0}/{1}.pvd".format(self.output_folder, file_name)) for file_name in self._wanted_files[func_group]},
        )

    def diffuse_real_scale_test(self, initial_condition, start_time, end_time, time_step, output_time):
        # Ensure we have a land attr set
        if not hasattr(self, "land"):
            raise ValueError("No land geometry function found")
        sl_time = fd.Constant(25 * fd.sin(start_time / 100000 / 180 * math.pi))

        # Initialise our functions and out files
        funcs, out_files = self._initialise_features("real_scale")
        funcs["slf"] = fd.Function(self.function_space, name="sea_level")
        funcs["sea_level"] = fd.Function(
            self.function_space,
            val=fd.interpolate(sl_time, self.function_space),
            name="sea_level",
        )
        # if we don't do this, the sea level ouput has a random function name.
        # Not helpful...
        funcs["slf"].interpolate(funcs["sea_level"])

        out_files["land"].write(self.land)

        # Assign our two sediment buildup functions to the initial condition
        funcs["sed"].assign(initial_condition)
        funcs["sed_old"].assign(initial_condition)

        # Perform initial interpolation
        funcs["surface"].interpolate((
            ((self.land + funcs["sed"]) + self.land)
            + abs((self.land + funcs["sed"]) - self.land)
        ) / 2)
        funcs["thickness"].interpolate(funcs["surface"] - self.land)
        funcs["limiter"].interpolate(
            (funcs["surface"] - self.land)
            / (funcs["surface"] - self.land + TINY)
        )
        funcs["depth"].interpolate(funcs["sea_level"] - funcs["surface"])
        funcs["diff"].project((
            2 / (fd.sqrt(2 * math.pi))
            * fd.exp(-0.5 * funcs["depth"] ** 2)
        ) + 0.2022)

        # Specify our function to solve (ie. the diffusion equation)
        F = (
            fd.inner(
                (funcs["sed"] - funcs["sed_old"]) / time_step,
                self.test_function
            ) + funcs["limiter"]
            * funcs["diff"]
            * fd.inner(
                fd.grad(funcs["sed"] + self.land),
                fd.grad(self.test_function)
            )
        ) * fd.dx

        # Write out initial data to our blank output files
        out_files["surfaces"].write(funcs["surface"], funcs["sed"], time=0)
        out_files["layer_data"].write(funcs["diff"], funcs["thickness"], funcs["depth"], time=0)
        out_files["sea_level"].write(funcs["slf"], time=0)

        # Begin looping until our end time has been reached
        current_time = start_time
        while current_time <= end_time:
            # Perform the solve
            fd.solve(F == 0, funcs["sed"])

            # Iterate
            funcs["sed_old"].assign(funcs["sed"])
            current_time += time_step

            # Advance the functions
            funcs["limiter"].interpolate(
                (funcs["surface"] - self.land)
                / (funcs["surface"] - self.land + TINY)
            )
            funcs["surface"].interpolate(
                (
                    ((self.land + funcs["sed"]) + self.land)
                    + abs((self.land + funcs["sed"]) - self.land)
                ) / 2
            )
            funcs["depth"].interpolate(funcs["sea_level"] - funcs["surface"])
            funcs["diff"].interpolate((
                2 / fd.sqrt(2 * math.pi)
                * fd.exp(-0.5 * ((funcs["depth"] - 5) / 10) ** 2)
            ) + 0.2022)
            funcs["thickness"].interpolate(funcs["surface"] - self.land)
            sl_time.assign(25 * fd.sin(current_time / 100000 * 180 / math.pi))
            funcs["sea_level"].interpolate(sl_time)
            funcs["slf"].interpolate(funcs["sea_level"])

            # If the current cycle is one where we should output, write to our files
            if current_time % output_time == 0:
                out_files["layer_data"].write(
                    funcs["diff"],
                    funcs["thickness"],
                    funcs["depth"],
                    time=current_time
                )
                out_files["surfaces"].write(funcs["surface"], funcs["sed"], time=current_time)
                out_files["sea_level"].write(funcs["slf"], time=current_time)


    def diffuse_carbonate_test(self, initial_condition, start_time, end_time, time_step, output_time):
        # Ensure we have a land attr set
        if not hasattr(self, "land"):
            raise ValueError("No land geometry function found")

        # Initialise our functions and out files
        sl_time = fd.Constant(25 * fd.sin(start_time / 1000000 / 180 * math.pi))
        funcs, out_files = self._initialise_features("carbonates")
        funcs["sea_level"] = fd.Function(
            self.function_space,
            val=fd.interpolate(sl_time, self.function_space),
            name="sea_level"
        )
        funcs["slf"].interpolate(funcs["sea_level"])  # Hack to prevent random file names

        # Assign our initial conditions
        funcs["sed"].assign(initial_condition)
        funcs["sed_old"].assign(initial_condition)
