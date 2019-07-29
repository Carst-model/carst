import copy
import math

import firedrake as fd

from .functions import FunctionContainer
from .functions import carst_funcs as f
from .options import CarstOptions
from .processes import (DIFFUSION_EQUATION_GENERIC, INIT_INTERPOLATION_ORDER,
                        advance_carbonates, advance_diffusion)


class CarstModel():
    """Simulates sediment formation on the seabed.
    """
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

    def __init__(self, options):
        """Initialise CarstModel instance.

        :param carst.options.CarstOptions options: The CarstOptions instance whose information you want to use to initialise the model.
        :raise: TypeError if options is not of type CarstOptions.
        :return: The initialised CarstModel instance.
        :rtype: carst.solver.CarstModel
        """
        # Type Checking
        if not isinstance(options, CarstOptions):
            raise TypeError(
                f"options must be of type CarstOptions, not {str(type(options))}"
            )

        # Get values from options
        self._options = options
        self._times = self._options["times"]
        self._out_files = self._options["out_files"]

        # Initialise function objects
        self._funcs = FunctionContainer(self._options, options["wanted_funcs"])

        # init sea level
        self._funcs[f.sea_level].assign(eval(self._options["sea_level"]))
        # Perform first output and interpolation
        self._funcs.interpolate(self._options, *INIT_INTERPOLATION_ORDER)

        # Initialise a diffusion equation if it's enabled, project diff_coeff
        if self._options["enabled_steps"]["diffusion"]:
            self._options["diffusion_equation"] = DIFFUSION_EQUATION_GENERIC(
                self._funcs, self._options)

        if self._options.get("initial_condition") is not None:
            self.set_condition(self._options["initial_condition"])

    @property
    def land(self):
        """:returns: The firedrake object containing the land function.
        :rtype: firedrake.function.Function
        """
        return self._options["land"]

    # Possibly change me to return by value/representation, not reference
    @property
    def funcs(self):
        """:returns: The FunctionContainer object which holds the current status of all the mathematical functions the model is using.
        :rtype: carst.functions.FunctionContainer
        """
        return self._funcs

    @property
    def coordinate_space(self):
        """:returns: The firedrake object describing the coordinate space the model is operating in.
        :rtype: firedrake.ufl.geometry.SpacialCoordinate
        """
        return self._options["coordinate_space"]

    @property
    def function_space(self):
        """:returns: The firedrake object describing the function space the model is operating in.
        :rtype: firedrake.functionspaceimpl.WithGeometry
        """
        return self._options["function_space"]

    @property
    def times(self):
        """:returns: A dict containing all the time variables the model is using.
        :rtype: dict
        """
        return copy.deepcopy(self._options["times"])

    @property
    def output_this_cycle(self):
        """:returns: True if model will write to output files this cycle.
        :rtype: bool
        """
        return self._times["current_time"] % self._times["output_time"] == 0

    def set_condition(self, condition):
        """Set the function describing the status of the sediment.

        :param firedrake.function.Function condition: the function to set the sediment to.
        """
        self._funcs[f.sed].assign(condition)
        self._funcs[f.sed_old].assign(condition)
        self._funcs.interpolate(self._options, *INIT_INTERPOLATION_ORDER)
        self._out_files.output(self._funcs, self._options)

    def advance(self):
        """Advance the simulation by a single time step.
        """
        # Increment time

        # Advance diffusion
        if self._options["enabled_steps"].get("diffusion"):
            advance_diffusion(self._funcs, self._options)

        # Advance carbonates
        if self._options["enabled_steps"].get("carbonates"):
            advance_carbonates(self._funcs, self._options)
            self._funcs[f.sed] += self._options['carbonate_production'] * self._funcs[f.light_attenuation]

        self._funcs[f.sed_old].assign(self._funcs[f.sed])

        # Output if necessary
        if self.output_this_cycle:
            print("At time step: " + str(self._times['current_time']))
            self._out_files.output(self._funcs, self._options)

        self._times["current_time"] += self._times["time_step"]

        # update sea level
        self._funcs[f.sea_level].assign(eval(self._options['sea_level']))
