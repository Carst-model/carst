import firedrake as fd
import enum
import math

TINY = 1e-10


# To add a function, add it to carst_funcs then add the corresponding logic to
# FunctionContainer._interpolation_funcs (make sure you get the key right!)
class carst_funcs(enum.Enum):
    sed = 1
    sed_old = 2
    surface = 3
    limiter = 4
    thickness = 5
    depth = 6
    diff = 7
    light_attenuation = 8
    sea_level = 9
    slf = 10


class FunctionContainer():
    _interpolation_funcs = {
        carst_funcs.surface: lambda land, funcs: (land + funcs[carst_funcs.sed] + land) + abs((land + funcs[carst_funcs.sed]) - land),
        carst_funcs.thickness: lambda land, funcs: (funcs[carst_funcs.surface] - land),
        carst_funcs.limiter: lambda land, funcs: (funcs[carst_funcs.surface] - land) / (funcs[carst_funcs.surface] - land + TINY),
        carst_funcs.depth: lambda land, funcs: funcs[carst_funcs.sea_level] - funcs[carst_funcs.surface],
        carst_funcs.diff: lambda land, funcs: 2 / fd.sqrt(2 * math.pi) * fd.exp(-0.5 * funcs[carst_funcs.depth] ** 2),
    }

    def __init__(self, solver, wanted_funcs):
        # Type checking
        for func in wanted_funcs:
            if not isinstance(func, carst_funcs):
                raise ValueError()

        self._solver = solver
        self._functions = {
            func_name: fd.Function(
                self._solver.function_space,
                name=str(func_name),
            ) for func_name in wanted_funcs
        }

    def __len__(self):
        return len(self._functions)

    def __getitem__(self, key):
        if not isinstance(key, carst_funcs):
            raise TypeError("The key must be a member of carst_funcs")
        return self._functions[key]

    def __iter__(self):
        return iter(self._functions.keys())

    def interpolate(self, function_names):
        if isinstance(function_names, carst_funcs):
            if function_names not in self.functions.keys():
                raise ValueError("That function isn't in this object")

            to_interpolate = [function_names]
        else:
            to_interpolate = function_names
            # Bait an exception if to_interpolate isn't iterable
            iter(to_interpolate)

        for function in to_interpolate:
            try:
                self.functions[function].interpolate(
                    FunctionContainer._interpolation_funcs[function](self._solver.land, self.functions)
                )
            except KeyError:
                pass
