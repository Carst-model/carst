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
    ),
    "diffusion": (
        f.sed,
        f.sed_old,
        f.limiter,
        f.diff_coeff,
    ),
    "carbonates": (f.light_attenuation, ),
}


def advance_diffusion(funcs, land, time_step, test_function) -> fd.Function:
    # "Copy" the sed function and solve
    sed_new = funcs[f.sed].copy(False)
    fd.solve(
        DIFFUSION_EQUATION_GENERIC(funcs, land, time_step, test_function) == 0,
        sed_new,
    )
    funcs.interpolate(*INTERPOLATION_ORDER)
    return sed_new


def advance_carbonates(funcs: FunctionContainer,
                       carbonate_production: float) -> fd.Function:
    funcs.interpolate(f.light_attenuation)
    return carbonate_production * funcs[f.light_attenuation]
