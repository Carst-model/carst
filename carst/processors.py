import firedrake as fd
from functions import carst_funcs

# Set numerical constants
TINY = 1e-10
DIFFUSION_EQUATION_GENERIC = lambda solver: (
    fd.inner(
        (solver.funcs[carst_funcs.sed] - solver.funcs[carst_funcs.sed_old]) / solver.time_step,
        solver.test_function,
    ) + (
        solver.funcs[carst_funcs.limiter]
        * solver.funcs[carst_funcs.diff_coeff]
        * fd.inner(
            fd.grad(solver.funcs[carst_funcs.sed] + solver.land),
            fd.grad(solver.test_function)
        )
    )
) * fd.dx

# Set interpolation order constants
INIT_INTERPOLATION_ORDER = (
    carst_funcs.surface,
    carst_funcs.thickness,
    carst_funcs.limiter,
    carst_funcs.depth,
    carst_funcs.sea_level,
)
INTERPOLATION_ORDER = (
    carst_funcs.limiter,
    carst_funcs.surface,
    carst_funcs.depth,
    carst_funcs.diff_coeff,
    carst_funcs.thickness,
)


def advance_diffusion(solver):
    # "Copy" the sed function and solve
    sed_new = solver.funcs[carst_funcs.sed].copy(True)
    fd.solve(
        DIFFUSION_EQUATION_GENERIC(solver) == 0,
        sed_new,
    )
    solver.funcs.interpolate(*INTERPOLATION_ORDER)
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
