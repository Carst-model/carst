from .functions import DIFF_COEFF_PROJECT, FunctionContainer
from .functions import carst_funcs as f
from .options import CarstOptions
from .processors import (INIT_INTERPOLATION_ORDER, advance_carbonates,
                         advance_diffusion)


class CarstModel():
    _WANTED_FILES = {
        "land":
        lambda solver: (solver.land),
        "layer_data":
        lambda solver: (
            solver.funcs[f.diff_coeff],
            solver.funcs[f.thickness],
            solver.funcs[f.depth],
        ),
        "surfaces":
        lambda solver: (
            solver.funcs[f.surface],
            solver.funcs[f.sed],
        ),
        "sea_level":
        lambda solver: (solver.funcs[f.sea_level]),
    }

    def __init__(self, options: CarstOptions):
        if not isinstance(options, CarstOptions):
            raise TypeError("Arg to CarstModel must be of type CarstOptions")

        useful_info = options.useful_data
        self.mesh, self.coordinate_space, self.function_space, self.test_function, self._sea_level_constant, self._land = useful_info[
            "workspace"]

        self.current_time, self._output_time, self._time_step = useful_info[
            "times"]

        self._funcs = FunctionContainer(self, useful_info["wanted_funcs"])
        print("Created the following functions:")
        print("\n".join((str(func) for func in self._funcs.items())))

        self.enabled_steps = useful_info["enabled_steps"]
        if self.enabled_steps["carbonates"]:
            self.carbonate_production = useful_info["optional"][0]
        self._out_files = useful_info["output_files"]

        self._out_files.output(self, ["land"])
        self.funcs.interpolate(*INIT_INTERPOLATION_ORDER)
        self._funcs[f.diff_coeff].project(DIFF_COEFF_PROJECT(self))

    @property
    def land(self):
        return self._land

    @property
    def funcs(self) -> FunctionContainer:
        return self._funcs

    @property
    def time_step(self) -> float:
        return self._time_step

    def set_condition(self, condition):
        self._funcs[f.sed].assign(condition)
        self._funcs[f.sed_old].assign(condition)

    def advance(self):
        self.current_time += self._time_step
        if self.enabled_steps.get("diffusion"):
            self._funcs[f.sed] = advance_diffusion(self)
        if self.enabled_steps.get("carbonates"):
            self._funcs[f.sed] = self._funcs[f.sed] + advance_carbonates(
                self._funcs, self.carbonate_production)

        if self.current_time % self._output_time == 0:
            self._out_files.output()
