import firedrake as fd
import enum
import math

TINY = 1e-10


class Available_Funcs(enum.Enum):
    sed = 1
    sed_old = 2
    surface = 3
    limiter = 4
    thickness = 5
    depth = 6
    diff = 7
    light_attenuation = 8


class Carst_Functions():
    _interpolation_funcs = {
        Available_Funcs.surface: lambda land, funcs: ((land + funcs["sed"]) + land) + abs((land + funcs["sed"]) - land),
        Available_Funcs.thickness: lambda land, funcs: (funcs["surface"] - land),
        Available_Funcs.limiter: lambda land, funcs: (funcs["surface"] - land) / (funcs["surface"] - land + TINY),
        Available_Funcs.depth: lambda land, funcs: funcs["sea_level"] - funcs["surface"],
        Available_Funcs.diff: lambda land, funcs: 2 / fd.sqrt(2 * math.pi) * fd.exp(-0.5 * funcs["depth"] ** 2),
    }

    def __init__(self, solver, wanted_funcs):
        if not set(wanted_funcs).issubset(set(Available_Funcs)):
            raise ValueError("Your wanted functions require an unsupported function")

        self._solver = solver
        self.functions = {
            fd.Function(
                self._solver.function_space,
                name=func_name
            ) for func_name in wanted_funcs
        }

    def interpolate(self, function_name):
        if function_name not in self.functions.keys():
            raise ValueError("That function isn't in this object")

        try:
            self.functions[function_name].interpolate(
                Carst_Functions._interpolation_funcs[function_name](self._solver.land, self.functions)
            )
        except KeyError:
            pass
