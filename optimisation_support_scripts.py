from carst_adjoint import *
from firedrake import Expression as firedrake_Expression # Use non-adjoint Expression

def select_distinct_gauges(det_xy, radius=100):
    # Choose a set of gauges less than a given radius from each other for optimisation
    n = len(det_xy)
    gauges_list = []
    # Assume gauges are appropriate
    data_filter = [True for i in range(n)]
    # Filter by distance to other gauges
    for i in range(n):
        j = i+1
        while j<n and data_filter[i]:
            if data_filter[j] and np.sqrt((det_xy[i][0] - det_xy[j][0])**2 + (det_xy[i][1] - det_xy[j][1])**2) < radius:
                # Points i and j are too close together: remove point i
                data_filter[i] = False
            j += 1

    gauges_list = [i for i, inc in enumerate(data_filter) if inc]
    return gauges_list

def generate_bump_functions(mesh2d, wells, gammas=None):
    P1_2d = FunctionSpace(mesh2d, 'CG', 1)  # Functionspace
    bump_functions = []
    i = 0
    # Default is 200m bump function size for all gauges
    if gammas is None:
        gammas = [200.0 for i in range(len(wells))]
    x = SpatialCoordinate(mesh2d)

    for xy in wells:
        x0 = xy[0]
        y0 = xy[1]
        gamma = gammas[i]
        this_bump_function = Function(P1_2d, name='bump_function_'+str(i))
        this_bump_function.interpolate(exp(-1/(gamma*gamma) * ((x[0]-x0)*(x[0]-x0) + (x[1]-y0)*(x[1]-y0))))
        this_bump_function.interpolate(this_bump_function/assemble(this_bump_function*dx))
        bump_functions.append(this_bump_function)
        i += 1

    return bump_functions
