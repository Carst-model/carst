# import firedrake as fd
from functions import carst_funcs
from processors import advance_carbonates, advance_diffusion
from options import CarstOptions


class CarstModel():
    _WANTED_FILES = {
        "land": lambda solver: (solver.land),
        "layer_data": lambda solver: (
            solver.funcs[carst_funcs.diff_coeff],
            solver.funcs[carst_funcs.thickness],
            solver.funcs[carst_funcs.depth],
        ),
        "surfaces": lambda solver: (
            solver.funcs[carst_funcs.surface],
            solver.funcs[carst_funcs.sed],
        ),
        "sea_level": lambda solver: (solver.funcs[carst_funcs.sea_level]),
    }

    def __init__(self, options):
        if not isinstance(options, CarstOptions):
            raise TypeError("Arg to CarstModel must be of type CarstOptions")

        useful_info = options.useful_info
        self.mesh, self.coordinate_space, self.function_space, self.test_function, self._sea_level_constant, self.land = useful_info[0]

        self.current_time, self._output_time, self._time_step = useful_info[1]

        self._funcs = useful_info[2]
        self.enabled_steps = useful_info[3]
        self._out_files = useful_info[4]

        self._out_files.output_selective("land")
        # Interpolate the funcs for the first time here!!!

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

    def advance(self):
        self.current_time += self._time_step
        if self.enabled_steps.get("diffusion"):
            self._funcs[carst_funcs.sed] = advance_diffusion(self)
        if self.enabled_steps.get("carbonates"):
            self._funcs[carst_funcs.sed] = self._funcs[carst_funcs.sed] + advance_carbonates(self._funcs, self.carbonate_production)

        if self.current_time % self._output_time == 0:
            self._out_files.output()
