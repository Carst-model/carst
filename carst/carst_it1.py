import math
import firedrake as fd

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
        self.land = fd.project(
            100 * fd.tanh(0.0005 * (self.coordinate_space[0] - 6000)),
            self.function_space,
            name="starting_topo",
        )

        # Declare the names of the functions we want for each solver
        self.wanted_functions = {
            "real_scale": (
                "sed",
                "sed_old",
                "surface",
                "limiter",
                "thickness",
                "depth",
                "diff",
            ),
        }

    def diffuse_real_scale_test(self, initial_condition, start_time, end_time, time_step, output_time):
        sl_time = fd.Constant(25 * fd.sin(start_time / 100000 / 180 * math.pi))

        # Initialise our dictionary of functions, declaring "slf" and "sea_level" on their own
        funcs = {func_name: fd.Function(self.function_space, name=func_name) for func_name in self.wanted_functions["real_scale"]}
        funcs["slf"] = fd.Function(self.function_space, name="sea_level")
        funcs["sea_level"] = fd.Function(
            self.function_space,
            val=fd.interpolate(sl_time, self.function_space),
            name="sea_level"
        )
        # if we don't do this, the sea level ouput has a random function name.
        # Not helpful...
        funcs["slf"].interpolate(funcs["sea_level"])

        outfile = fd.File(self.output_folder + "/land.pvd")
        outfile.write(self.land)

        # Assign our two sediment buildup functions to the initial condition
        funcs["sed"].assign(initial_condition)
        funcs["sed_old"].assign(initial_condition)

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
        wanted_files = ("surfaces", "layer_data", "sea_level")
        out_files = {file_name: fd.File("{0}/{1}.pvd".format(self.output_folder, file_name)) for file_name in wanted_files}
        out_files["surfaces"].write(funcs["surface"], funcs["sed"], time=0)
        out_files["layer_data"].write(funcs["diff"], funcs["thickness"], funcs["depth"], time=0)
        out_files["sea_level"].write(funcs["slf"], time=0)

        # Begin looping until our end time has been reached
        current_time = start_time
        while current_time <= end_time:
            fd.solve(F == 0, funcs["sed"])
            funcs["sed_old"].assign(funcs["sed"])
            current_time += time_step
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
                out_files["surface"].write(funcs["surface"], funcs["sed"], time=current_time)
                out_files["sea_level"].write(funcs["slf"], time=current_time)


# Run with test values
if __name__ == "__main__":
    my_solver = Diffuse_Solver(fd.RectangleMesh(50, 25, 10000, 5000), "output")
    my_solver.diffuse_real_scale_test(
        fd.project(
            (
                20000
                * (1 / (2 * fd.sqrt(2*math.pi*250*250)))
                * fd.exp(
                    -((my_solver.coordinate_space[0]-6000) * (my_solver.coordinate_space[0]-6000))
                    / (2 * 250 * 250)
                )
            )
            + (
                50000
                * (1 / (2 * fd.sqrt(2 * math.pi * 1000 * 1000)))
                * fd.exp(
                    -((my_solver.coordinate_space[0]-4000) * (my_solver.coordinate_space[0]-4000))
                    / (2 * 1000 * 1000)
                )
            ), my_solver.function_space),
        0,
        20000,
        50,
        500,
    )
