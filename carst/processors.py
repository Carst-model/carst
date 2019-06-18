import firedrake as fd
from functions import carst_funcs

# Set numerical constants
TINY = 1e-10
DIFFUSION_EQUATION_GENERIC = lambda solver: (
    fd.inner(
        (solver.funcs[carst_funcs.sed] - solver.funcs[carst_funcs.sed_old]) / solver.time_step,
        solver.test_function
    ) + solver.funcs[carst_funcs.limiter]
    * solver.funcs[carst_funcs.diff_coeff]
    * fd.inner(
        fd.grad(solver.funcs[carst_funcs.sed] + solver.land),
        fd.grad(solver.test_function)
    )
) * fd.dx


def advance_diffusion(solver):
    sed_new = fd.Function(solver.function_space, name="sed")
    fd.solve(
        DIFFUSION_EQUATION_GENERIC(solver) == 0,
        sed_new,
    )
    return sed_new


def advance_carbonates(funcs, carbonate_production):
    funcs.interpolate(carst_funcs.light_attenuation)
    return carbonate_production * funcs[carst_funcs.light_attenuation]


PROCESSOR_NEEDED_FUNCS = {
    "diffusion": (
        carst_funcs.sed,
        carst_funcs.sed_old,
        carst_funcs.limiter,
        carst_funcs.diff_coeff,
    ),
    "carbonates": (
        carst_funcs.light_attenuation,
    ),
}
