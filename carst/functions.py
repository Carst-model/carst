import firedrake as fd
import math

TINY = 1e-10


class Carst_Functions():
    available_funcs = [
        "sed",
        "sed_old",
        "surface",
        "limiter",
        "thickness",
        "depth",
        "diff",
        "light_attenuation",
    ]

    _interpolation_funcs = {
        "surface": lambda land, funcs: ((land + funcs["sed"]) + land) + abs((land + funcs["sed"]) - land),
        "thickness": lambda land, funcs: (funcs["surface"] - land),
        "limiter": lambda land, funcs: (funcs["surface"] - land) / (funcs["surface"] - land + TINY),
        "depth": lambda land, funcs: funcs["sea_level"] - funcs["surface"],
        "diff": lambda land, funcs: 2 / fd.sqrt(2 * math.pi) * fd.exp(-0.5 * funcs["depth"] ** 2),
    }

    def __init__(self, solver, wanted_funcs):
        if not set(wanted_funcs).issubset(set(Carst_Functions.available_funcs)):
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

        self.functions[function_name].interpolate(
            Carst_Functions._interpolation_funcs[function_name](self._solver.land, self.functions)
        )
