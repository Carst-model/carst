import firedrake as fd
from functions import carst_funcs, FunctionContainer
from processors import PROCESSOR_NEEDED_FUNCS, advance_carbonates, advance_diffusion


class CarstModel():
    _WANTED_FILES = {
        "land": lambda land: (land),
        "layer_data": lambda funcs: (funcs[carst_funcs.diff_coeff], funcs[carst_funcs.thickness], funcs[carst_funcs.depth]),
        "surfaces": lambda funcs: (funcs[carst_funcs.surface], funcs[carst_funcs.sed]),
        "sea_level": lambda funcs: (funcs[carst_funcs.sea_level]),
    }

    def __init__(self, base_mesh, land, sea_level_constant, times, **kw_args):
        if not isinstance(base_mesh, fd.mesh.MeshGeometry):
            raise TypeError("base_mesh not of type firedrake.Mesh")
        if not isinstance(sea_level_constant, fd.Constant):
            raise TypeError("sea_level_constant not of type firedrake.Constant")

        # Need to check whether the output folder actually exists with os

        # Store the passed values
        self._sea_level_constant = sea_level_constant
        self.current_time, self._time_step, self._output_time = times
        self._output_folder = kw_args.get("output_folder")
        self.mesh = base_mesh

        # Mark the steps in the process we want
        self.enabled_steps = {
            "diffusion": kw_args.get("diffusion"),
            "carbonates": kw_args.get("carbonates"),
        }

        # Generate our workspace from the mesh
        self.coordinate_space = fd.SpatialCoordinate(self.mesh)
        self.function_space = fd.FunctionSpace(self.mesh, "CG", 1)
        self.test_function = fd.TestFunction(self.function_space)

        # Get our land, based in our workspace
        self._land = land(self.coordinate_space, self.function_space)

        # Initialise the funcs we need
        self._funcs = FunctionContainer(
            self,
            {
                func_name for processor in PROCESSOR_NEEDED_FUNCS.values() for func_name in processor
            },
        )
        # self._funcs.interpolate(*type(self)._INIT_INTERPOLATION_ORDER)

        # Initialise _out_files and write land to them (if we have an output passed)
        if self._output_folder is not None:
            self._out_files = {
                file_name: fd.File("{0}/{1}.pvd".format(self._output_folder, file_name)) for file_name in type(self)._WANTED_FILES.keys()
            }
            self._out_files["land"].write(self._land)

    @property
    def land(self):
        return self._land

    @property
    def funcs(self):
        return self._funcs

    @property
    def time_step(self):
        return self._time_step

    def set_condition(self, condition):
        self._funcs[carst_funcs.sed].assign(condition)
        self._funcs[carst_funcs.sed_old].assign(condition)

    def output(self):
        for file_name in type(self)._WANTED_FILES:
            if file_name != "land":
                self._out_files[file_name].write(
                    *type(self)._WANTED_FILES[file_name](self._funcs),
                    time=self.current_time,
                )

    def advance(self):
        self.current_time += self._time_step
        if self.enabled_steps.get("diffusion"):
            self._funcs[carst_funcs.sed] = advance_diffusion(self)
        if self.enabled_steps.get("carbonates"):
            self._funcs[carst_funcs.sed] = self._funcs[carst_funcs.sed] + advance_carbonates(self._funcs, self.carbonate_production)

        if self.current_time % self._output_time == 0:
            self.output()
