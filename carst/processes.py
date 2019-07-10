import firedrake as fd

from .functions import FunctionContainer
from .functions import carst_funcs as f

# Set numerical constants
TINY = 1e-10


def DIFFUSION_EQUATION_GENERIC(funcs: FunctionContainer,
                               options) -> fd.Function:
    return (fd.inner(
        (funcs[f.sed] - funcs[f.sed_old]) / options["times"]["time_step"],
        options["test_function"],
    ) + funcs[f.limiter] * options['diff_coeff'] * #funcs[f.diff_coeff] *
            fd.inner(fd.grad(funcs[f.sed] + options["land"]),
                     fd.grad(options["test_function"]))) * fd.dx


# Set interpolation order constants
INIT_INTERPOLATION_ORDER = (
    #f.sea_level,
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
    fd.solve(
        options['diffusion_equation'] == 0,
        funcs[f.sed]
    )
    funcs[f.sed_old].assign(funcs[f.sed])
    funcs.interpolate(options, *INTERPOLATION_ORDER)


def advance_carbonates(funcs: FunctionContainer, options) -> fd.Function:
    funcs.interpolate(options, f.light_attenuation)
    return options["carbonate_production"] * funcs[f.light_attenuation]
