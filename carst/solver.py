import copy

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
        self._options = options
        self._times = self._options["times"]

        self._funcs = FunctionContainer(self, options["wanted_funcs"])
        print("Created the following functions:")
        print("\n".join((str(func) for func in self._funcs.items())))

        self._out_files = self._options["out_files"]

        self._out_files.output(self._funcs, self._options["land"],
                               self._times["current_time"], ("land", ))
        self._funcs.interpolate(self._options["land"],
                                self._options["sea_level_constant"],
                                *INIT_INTERPOLATION_ORDER)
        self._funcs[f.diff_coeff].project(DIFF_COEFF_PROJECT(self._funcs))

    @property
    def land(self):
        return self._options["land"]

    @property
    def funcs(self) -> FunctionContainer:
        return self._funcs

    @property
    def coordinate_space(self):
        return self._options["coordinate_space"]

    @property
    def function_space(self):
        return self._options["function_space"]

    @property
    def times(self) -> float:
        return copy.deepcopy(self._options["times"])

    def set_condition(self, condition):
        self._funcs[f.sed].assign(condition)
        self._funcs[f.sed_old].assign(condition)

    def advance(self):
        self._times["current_time"] += self._times["time_step"]
        if self._options["enabled_steps"].get("diffusion"):
            self._funcs[f.sed] = advance_diffusion(
                self._funcs, self._options["land"], self._times["time_step"],
                self._options["test_function"])
        if self._options["enabled_steps"].get("carbonates"):
            self._funcs[f.sed] = self._funcs[f.sed] + advance_carbonates(
                self._funcs, self._options["carbonate_production"])

        if self._times["current_time"] % self._times["output_time"] == 0:
            self._out_files.output()
