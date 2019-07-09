import enum
import math
from collections import UserDict
from typing import Iterable

import firedrake as fd

TINY = 1e-10


# To add a function, add it to carst_funcs then add the corresponding logic to
# FunctionContainer._INTERPOLATION_FUNCS (make sure you get the key right!)
class carst_funcs(enum.Enum):
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
    _INTERPOLATION_FUNCS = {
        carst_funcs.surface:
        lambda land, funcs, sea_level_constant: ((land + funcs[
            carst_funcs.sed] + land) + abs((land + funcs[carst_funcs.sed]) -
                                           land)),
        carst_funcs.thickness:
        lambda land, funcs, sea_level_constant: (funcs[carst_funcs.surface] -
                                                 land),
        carst_funcs.limiter:
        lambda land, funcs, sea_level_constant: ((funcs[
            carst_funcs.surface] - land) / (funcs[carst_funcs.surface] - land +
                                            TINY)),
        carst_funcs.depth:
        lambda land, funcs, sea_level_constant: (funcs[carst_funcs.sea_level] -
                                                 funcs[carst_funcs.surface]),
        carst_funcs.diff_coeff:
        lambda land, funcs, sea_level_constant: (2 / fd.sqrt(
            2 * math.pi) * fd.exp(-0.5 * funcs[carst_funcs.depth]**2)),
        carst_funcs.sea_level:
        lambda land, funcs, sea_level_constant: sea_level_constant,
    }

    def __init__(self, solver, wanted_funcs: Iterable[carst_funcs]):
        function_space = solver.function_space
        super().__init__({
            func_name: fd.Function(
                function_space,
                name=str(func_name),
            )
            for func_name in wanted_funcs
        })

    def __getitem__(self, key: carst_funcs) -> fd.Function:
        if not isinstance(key, carst_funcs):
            raise TypeError("Key not a member of carst_funcs")
        return super().__getitem__(key)

    def __setitem__(self, key: carst_funcs, val: fd.Function):
        if not isinstance(key, carst_funcs):
            raise TypeError("Key not a member of carst_funcs")
        if not isinstance(val, fd.Function):
            raise TypeError("Value not of type firedrake.Function")
        super().__setitem__(key, val)

    def interpolate(
            self,
            options,
            *function_names: Iterable[carst_funcs],
    ):
        for name in function_names:
            try:
                self[name].interpolate(
                    FunctionContainer._INTERPOLATION_FUNCS[name](
                        options["land"], self, options["sea_level_constant"]))
            except KeyError:
                continue


def DIFF_COEFF_PROJECT(funcs: FunctionContainer) -> fd.Function:
    return (((2 / fd.sqrt(2 * math.pi)) *
             fd.exp(-0.5 * funcs[carst_funcs.depth]**2)) + 0.2022)
