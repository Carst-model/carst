import firedrake as fd

from .functions import FunctionContainer
from .functions import carst_funcs as f

# Set numerical constants
TINY = 1e-10


def DIFFUSION_EQUATION_GENERIC(funcs: FunctionContainer, land,
                               time_step: float, test_function) -> fd.Function:
    return (fd.inner(
        (funcs[f.sed] - funcs[f.sed_old]) / time_step,
        test_function,
    ) + funcs[f.limiter] * funcs[f.diff_coeff] * fd.inner(
        fd.grad(funcs[f.sed] + land), fd.grad(test_function))) * fd.dx


# Set interpolation order constants
INIT_INTERPOLATION_ORDER = (
    f.sea_level,
    f.surface,
    f.thickness,
    f.limiter,
    f.depth,
)
INTERPOLATION_ORDER = (
    f.limiter,
    f.surface,
    f.depth,
    f.diff_coeff,
    f.thickness,
)

PROCESSOR_NEEDED_FUNCS = {
    "basic": (
        f.sed,
        f.sed_old,
        f.surface,
        f.thickness,
        f.depth,
        f.sea_level,
    ),
    "diffusion": (
        f.sed,
        f.sed_old,
        f.limiter,
        f.diff_coeff,
    ),
    "carbonates": (f.light_attenuation, ),
}


def advance_diffusion(funcs: FunctionContainer, options):
    # "Copy" the sed function and solve
    fd.solve(
        DIFFUSION_EQUATION_GENERIC(funcs, options["land"],
                                   options["times"]["time_step"],
                                   options["test_function"]) == 0,
        funcs[f.sed],
    )
    funcs.interpolate(options, *INTERPOLATION_ORDER)
    funcs[f.sed_old] = funcs[f.sed]


def advance_carbonates(funcs: FunctionContainer, options) -> fd.Function:
    funcs.interpolate(options, (f.light_attenuation, ))
    return options["carbonate_production"] * funcs[f.light_attenuation]
