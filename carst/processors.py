import firedrake as fd
from functions import carst_funcs as f

# Set numerical constants
TINY = 1e-10
DIFFUSION_EQUATION_GENERIC = lambda solver: (fd.inner(
    (solver.funcs[f.sed] - solver.funcs[f.sed_old]) / solver.time_step,
    solver.test_function,
) + (solver.funcs[f.limiter] * solver.funcs[f.diff_coeff] * fd.inner(
    fd.grad(solver.funcs[f.sed] + solver.land), fd.grad(solver.test_function)))
                                             ) * fd.dx

# Set interpolation order constants
INIT_INTERPOLATION_ORDER = (
    f.surface,
    f.thickness,
    f.limiter,
    f.depth,
    f.sea_level,
)
INTERPOLATION_ORDER = (
    f.limiter,
    f.surface,
    f.depth,
    f.diff_coeff,
    f.thickness,
)


def advance_diffusion(solver):
    # "Copy" the sed function and solve
    sed_new = solver.funcs[f.sed].copy(True)
    fd.solve(
        DIFFUSION_EQUATION_GENERIC(solver) == 0,
        sed_new,
    )
    solver.funcs.interpolate(*INTERPOLATION_ORDER)
    return sed_new


def advance_carbonates(funcs, carbonate_production):
    funcs.interpolate(f.light_attenuation)
    return carbonate_production * funcs[f.light_attenuation]


PROCESSOR_NEEDED_FUNCS = {
    "diffusion": (
        f.sed,
        f.sed_old,
        f.limiter,
        f.diff_coeff,
    ),
    "carbonates": (f.light_attenuation, ),
}
