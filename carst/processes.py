from .firedrake import *

from .functions import FunctionContainer
from .functions import carst_funcs as f

# Set numerical constants
TINY = 1e-10


def DIFFUSION_EQUATION_GENERIC(funcs, options):
    """Generate a firedrake object representing the diffusion equation for the current model.

    :param carst.functions.FunctionContainer funcs: The *FunctionContainer* instance the model is working from.
    :param carst.options.CarstOptions options: The *CarstOptions* instance the model is working from.
    :returns: A firedrake Function to be used for solving the diffusion equation with the current model (reusable).
    :rtype: firedrake.function.Function
    """
    return (inner(
        (funcs[f.sed] - funcs[f.sed_old]) / options["times"]["time_step"],
        options["test_function"],
    ) + funcs[f.limiter] * options['diff_coeff'] *  #funcs[f.diff_coeff] *
            inner(grad(funcs[f.sed] + options["land"]),
                     grad(options["test_function"]))) * dx


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


def advance_diffusion(funcs, options):
    """Perform diffusion simulation over 1 time step.

    Note that this also performs interpolation over the rest of the *FunctionContainer* since it is the most basic function of the model.

    :param carst.functions.FunctionContainer funcs: The *FunctionContainer* instance the model is working from. This is modified in-place.
    :param carst.options.CarstOptions options: The *CarstOptions* instance the model is working from.
    """
    solve(options['diffusion_equation'] == 0, funcs[f.sed])
    funcs[f.sed_old].assign(funcs[f.sed])
    funcs.interpolate(options, *INTERPOLATION_ORDER)


def advance_carbonates(funcs, options):
    """Perform carbonate simulation over 1 time step.

    :param carst.functions.FunctionContainer funcs: The *FunctionContainer* instance the model is working from. This is modified in-place.
    :param carst.options.CarstOptions options: The *CarstOptions* instance the model is working from.
    """
    funcs.interpolate(options, f.light_attenuation)
    return options["carbonate_production"] * funcs[f.light_attenuation]
