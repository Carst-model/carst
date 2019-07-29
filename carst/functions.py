import enum
import math
from collections import UserDict
from typing import Sequence

import firedrake as fd


# To add a function, add it's label to carst_funcs then add the corresponding logic to
# FunctionContainer._INTERPOLATION_FUNCS (make sure you get the key right!)
class carst_funcs(enum.Enum):
    """Enumerator to correspond to the functions the model relies on.
    """
    sed = 1
    sed_old = 2
    surface = 3
    limiter = 4
    thickness = 5
    depth = 6
    diff_coeff = 7
    light_attenuation = 8
    sea_level = 9


class FunctionContainer(UserDict):
    """Mapping to hold functions required for modelling.

    Inherits *dict*/*UserDict* to allow easy access to parent data. Note that only `carst.functions.carst_funcs`_ values are acceptable as keys.

    :param carst.options.CarstOptions options: The CarstOptions instance the model is working off. Needed for supplementary data.
    :param Sequence[carst.functions.carst_funcs] wanted_funcs: The collection of carst_funcs which specifies which function types the object should contain.
    :returns: An initialised carst.functions.FunctionContainer instance. The functions are still blank (ie. not interpolated).
    """
    # Workaround for python's lack of a switch/match statement >:(
    _INTERPOLATION_FUNCS = {
        carst_funcs.surface:
        lambda funcs, options: ((((2.0 * options["land"]) + funcs[
            carst_funcs.sed]) + abs(funcs[carst_funcs.sed])) / 2.0),
        carst_funcs.thickness:
        lambda funcs, options: (funcs[carst_funcs.surface] - options["land"]),
        carst_funcs.limiter:
        lambda funcs, options:
        ((funcs[carst_funcs.surface] - options["land"]) /
         (funcs[carst_funcs.surface] - options["land"] + 1e-10)),
        carst_funcs.depth:
        lambda funcs, options:
        (funcs[carst_funcs.sea_level] - funcs[carst_funcs.surface]),
        carst_funcs.diff_coeff:
        lambda funcs, options: (options["diff_coeff"] * ((2. / fd.sqrt(2. * math.pi)) * fd.exp(-0.5 * ((funcs[
            carst_funcs.depth]-5.0)/10.0)**2))),
        carst_funcs.light_attenuation:
        lambda funcs, options: (1.0/(1.+fd.exp(-2*funcs[carst_funcs.depth]*25.0))* fd.exp(-1.0*funcs[carst_funcs.depth]/10.)),

        #carst_funcs.sea_level:
        #lambda funcs, options: options["sea_level"],
    }

    def __init__(self, options, wanted_funcs: Sequence[carst_funcs]):
        function_space = options["function_space"]
        super().__init__({
            func_name: fd.Function(
                function_space,
                name=str(func_name),
            )
            for func_name in wanted_funcs
        })

    # Enforce type checking on __getitem__ and __setitem__
    def __getitem__(self, key):
        if not isinstance(key, carst_funcs):
            raise TypeError(f"Key {str(key)} not a member of carst_funcs")
        return super().__getitem__(key)

    def __setitem__(self, key, val):
        if not isinstance(key, carst_funcs):
            raise TypeError(f"Key {str(key)} not a member of carst_funcs")
        if not isinstance(val, fd.Function):
            raise TypeError(f"Value {str(val)} not of type firedrake.Function")
        super().__setitem__(key, val)

    def interpolate(self, options, *function_names):
        """Interpolate the functions contained in the module by one time step.

        :param carst.options.CarstOptions options: The options set that the model is currently working on.
        :param Sequence[carst.functions.carst_funcs] function_names: The functions to be interpolated.
        """
        for name in function_names:
            try:
                self[name].interpolate(
                    FunctionContainer._INTERPOLATION_FUNCS[name](self,
                                                                 options))
            except KeyError:
                continue

    def __repr__(self):
        return str(__class__).split("'")[1] + "(" + ", ".join(
            self.keys()) + ")"
