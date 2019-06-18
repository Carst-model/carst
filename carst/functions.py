import firedrake as fd
import enum
import math

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


class FunctionContainer():
    _INTERPOLATION_FUNCS = {
        carst_funcs.surface: lambda solver, funcs: (solver.land + funcs[carst_funcs.sed] + solver.land) + abs((solver.land + funcs[carst_funcs.sed]) - solver.land),
        carst_funcs.thickness: lambda solver, funcs: (funcs[carst_funcs.surface] - solver.land),
        carst_funcs.limiter: lambda solver, funcs: (funcs[carst_funcs.surface] - solver.land) / (funcs[carst_funcs.surface] - solver.land + TINY),
        carst_funcs.depth: lambda solver, funcs: funcs[carst_funcs.sea_level] - funcs[carst_funcs.surface],
        carst_funcs.diff_coeff: lambda solver, funcs: 2 / fd.sqrt(2 * math.pi) * fd.exp(-0.5 * funcs[carst_funcs.depth] ** 2),
        carst_funcs.sea_level: lambda solver, funcs: solver.sea_level_constant,
    }

    def __init__(self, solver, wanted_funcs):
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
            raise TypeError("Key not a member of carst_funcs")
        return self._functions[key]

    def __setitem__(self, key, val):
        if not isinstance(key, carst_funcs):
            raise TypeError("Key not a member of carst_funcs")
        if not isinstance(val, fd.Function):
            raise TypeError("Value not of type firedrake.Function")
        self._functions[key] = val

    def __iter__(self):
        return iter(self._functions.keys())

    def add(self, function_name, function):
        if not isinstance(function_name, carst_funcs):
            raise TypeError("function_name must be a carst_funcs instance")
        if not isinstance(function, fd.Function):
            raise TypeError("function must be of type firedrake.Function")
        if function_name in self._functions.keys():
            raise IndexError("There is already a copy of that function")
        self._functions[function_name] = function

    def interpolate(self, *function_names):
        for name in function_names:
            try:
                self._functions[name].interpolate(
                    FunctionContainer._INTERPOLATION_FUNCS[name](self._solver, self._functions)
                )
            except KeyError:
                pass
