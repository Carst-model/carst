import copy
import firedrake as fd
from .functions import FunctionContainer
from .functions import carst_funcs as f
from .options import CarstOptions
from .processes import (DIFFUSION_EQUATION_GENERIC, INIT_INTERPOLATION_ORDER,
                        advance_carbonates, advance_diffusion)


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
        # Type Checking
        if not isinstance(options, CarstOptions):
            raise TypeError("Arg to CarstModel must be of type CarstOptions")

        # Get values from options
        self._options = options
        self._times = self._options["times"]
        self._out_files = self._options["out_files"]

        # Initialise function objects
        self._funcs = FunctionContainer(self, options["wanted_funcs"])

        # Perform first output and interpolation
        self._funcs.interpolate(self._options, *INIT_INTERPOLATION_ORDER)

        # Initialise a diffusion equation if it's enabled, project diff_coeff
        if self._options["enabled_steps"]["diffusion"]:
            self._options["diffusion_equation"] = DIFFUSION_EQUATION_GENERIC(
                self._funcs, self._options)

    @property
    def land(self):
        return self._options["land"]

    # Possibly change me to return by value/representation, not reference
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
    def times(self) -> dict:
        return copy.deepcopy(self._options["times"])

    def set_condition(self, condition):
        self._funcs[f.sed].assign(condition)
        self._funcs[f.sed_old].assign(condition)
        self._funcs.interpolate(self._options, *INIT_INTERPOLATION_ORDER)        
        self._out_files.output(self._funcs, self._options)


    @property
    def output_this_cycle(self):
        return self._times["current_time"] % self._times["output_time"] == 0

    def advance(self):
        # Increment time

        # Advance diffusion
        if self._options["enabled_steps"].get("diffusion"):
            advance_diffusion(self._funcs, self._options)

        # Advance carbonates
        if self._options["enabled_steps"].get("carbonates"):
            self._funcs[f.sed] = self._funcs[f.sed] + advance_carbonates(
                self._funcs, self._options)

        # Output if necessary
        if self.output_this_cycle:
            print("At time step: "+str(self._times['current_time']))
            self._out_files.output(self._funcs, self._options)

        self._times["current_time"] += self._times["time_step"]

        # update sea level
        self._funcs[f.sea_level].assign(25*fd.sin(self._times['current_time']/100000*180/3.14159))
