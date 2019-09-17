#!/usr/bin/env python3
import copy
import math
import random
import numpy as np
import h5py
from pyadjoint.optimization.optimization import minimise
import optimisation_support_scripts
from firedrake_adjoint import *




START_TIME = 0
OUTPUT_TIME = 500
TIME_STEP = 50
OUTPUT_FOLDER = "output"

angle = math.radians(1)
run_time = 50000.

mesh = RectangleMesh(200, 200, 10000, 10000)
# set up forward model

# Example land function generator
def EXAMPLE_LAND(coordinate_space, function_space):
    return project(100 * tanh(0.0005 * (coordinate_space[0] - 6000) ),
                      function_space,
                      name="starting_topo")


# Example initial condition
def EXAMPLE_INITIAL_COND(coordinate_space, function_space):
    return project(
        (20000 * (1 / (2 * sqrt(2 * math.pi * 250 * 250))) * exp(-(
            (coordinate_space[0] - 6000) *
            (coordinate_space[0] - 6000)) / (2 * 250 * 250))) +
        (50000 * (1 / (2 * sqrt(2 * math.pi * 1000 * 1000))) * exp(-(
            (coordinate_space[0] - 4000) *
            (coordinate_space[0] - 4000)) / (2 * 1000 * 1000))),
        function_space)



# Initialise a solver and add land
my_options = CarstOptions(
    initialisation_method.raw_values,
    mesh,
    EXAMPLE_LAND,
    "25 * sin(t / 50000 * 180 / 3.142)",
    (
        START_TIME,
        TIME_STEP,
        OUTPUT_TIME,
        run_time
    ),
    output_folder=OUTPUT_FOLDER,
    diffusion=True,
    carbonates=True,
    diff_coeff=1.0,
    carbonate_production = 3.0,
)
my_solver_real_scale = CarstModel(my_options)

# DiffuseSolver doesn't play nice with deepcopy apparently
# my_solver_carbonates = copy.deepcopy(my_solver_real_scale)

# Run with a sample initial condition
my_solver_real_scale.set_condition(
    EXAMPLE_INITIAL_COND(my_solver_real_scale.coordinate_space,
                         my_solver_real_scale.function_space))

# add wells
# random locations between (0,0) and (6000,10000), which equates to where sed gets deposited
# also specify a minimum distance of 100m (~2 cells)
#coords = np.random.uniform(0.0, 6000., size=(25,2))
#DEBUG
coords = np.array([[4284.96714447, 4171.96162372],
 [5474.81825731,  957.02820747],
 [5270.17849967, 4238.39867382],
 [4085.80085666, 4083.45236598],
 [2221.11864083, 4613.83299363],
 [4644.53683661, 3753.25072896],
 [1305.47062266, 1793.15854649],
 [5401.02092879, 3953.27350635],
 [ 442.9601122,  4072.01711454],
 [3532.5004372,   585.54564048],
 [ 360.32187667, 2010.46646992],
 [2673.93267268, 2588.49886704],
 [ 749.00488793, 4034.77676085],
 [1076.52347337, 5278.88237022],
 [5373.76125456,  488.16742656],
 [ 738.79240511, 2649.15870396],
 [5414.72854895,  823.56587848],
 [2412.55908675, 5171.62168942],
 [2207.27667454,   29.42633511],
 [1624.22805545, 3404.75201571],
 [3852.06939361, 3387.36797005],
 [3311.00284404,  321.09196699],
 [5824.60923734, 2816.45891231],
 [4198.22700678, 3861.86727105],
 [4469.80783981, 4645.32027827]])
# a fixed set of coords whilst I debug. These are random though!

#coords = np.array([[0,0],[1000,1000]])

cb = callback.DetectorsCallback(my_solver_real_scale, 
                                      coords, [f.sed] ,"palaeosed", outputdir="output")
my_solver_real_scale.callbacks.add(cb, 'timestep')

# Iterate
#my_solver_real_scale.iterate()


# we now have the wells in the h5 detector file in outputs and these contain the palaeo depth at each timestep

# now we construct the adjoint
# We'll create a new solver object and intialise with an incorrect, spatially constant carb production. 
# but the optimisation can vary this value spatially.
# We'll then optimise and hopefully find a value of 3.0 across the whole domain
# Initialise a solver and add land
diffusion = Function(FunctionSpace(mesh, "CG", 1), name='diffusion')
diffusion.assign(5.0)
#carb_production = Constant(0.0)

my_options2 = CarstOptions(
    initialisation_method.raw_values,
    mesh,
    EXAMPLE_LAND,
    "25 * sin(t / 50000 * 180 / 3.142)",
    (
        START_TIME,
        TIME_STEP,
        OUTPUT_TIME,
        run_time
    ),
    output_folder="output2",
    diffusion=True,
    carbonates=True,
    diff_coeff=diffusion,
    carbonate_production = 3.0,
)
solver_obj = CarstModel(my_options2)
solver_obj.set_condition(
    EXAMPLE_INITIAL_COND(solver_obj.coordinate_space,
                         solver_obj.function_space))
# these are our well from the initial forward run
ff = h5py.File("output/diagnostic_palaeodepth.hdf5", 'r')
wells_from_file = {}
global times

times = np.array(ff['time'])
times = times.reshape(times.shape[0])
ndets = len(coords)
fill = len(str(ndets))
for n in range(ndets):
    detector_name = 'detector{:0{fill}d}'.format(n, fill=fill)
    wells_from_file[n] = np.array(ff[detector_name])


# a place to store our palaeodepth at the end of the simulation
pd_gauge_fields = [Constant(0.0) for i in range(len(coords))]
# Create bump functions and output to files. Our wells don't align exactly to the mesh
bump_functions = optimisation_support_scripts.generate_bump_functions(mesh, coords)
[File('output2/bump_'+str(i)+'.pvd').write(b) for i, b in enumerate(bump_functions)]
# create my functional
functional_expression = 0
for i in range(len(coords)):
    functional_expression += inner(solver_obj._funcs[f.sed]-pd_gauge_fields[i], solver_obj._funcs[f.sed]-pd_gauge_fields[i]) * bump_functions[i] * dx

# Initialise functional to zero
global time_integrated_functional
time_integrated_functional = 0.0


def update_forcings(t):
    # this is where we're going to update things each time step

    # perform some adjoint magic
    bv = solver_obj._funcs[f.sed].create_block_variable()
    bv.checkpoint = solver_obj._funcs[f.sed]._ad_create_checkpoint()

    global times
    if t == 0:
        index = 0
    else:
        index = np.argmax(np.where(times <= t))
    sum_diff = 0
    # update well_from_file
    for i in range(len(coords)):
        pd_gauge_fields[i].assign(wells_from_file[i][index][0])

    global time_integrated_functional
    time_integrated_functional += assemble(functional_expression)*TIME_STEP
    


# run as normal
solver_obj.iterate(update_forcings=update_forcings)


# now we set up the adjoint
c = Control(diffusion)
djd0 = compute_gradient(time_integrated_functional, c)
print("Gradient = ",float(djd0))

rf = ReducedFunctional(time_integrated_functional, c)



td_opt = minimise(rf, bounds=[0, 10.0],
                  options={'maxiter': 50, 'pgtol': 1e-3})
File('optimal_carb_prod.pvd').write(td_opt)
