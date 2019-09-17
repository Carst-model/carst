#!/usr/bin/env python3
import copy
import math
import random
import numpy as np
import h5py
from pyadjoint.optimization.optimization import minimise
import optimisation_support_scripts
from firedrake import *
from firedrake_adjoint import *




START_TIME = 0
OUTPUT_TIME = 500
timestep = 50
OUTPUT_FOLDER = "output"

angle = math.radians(1)
run_time = 50000.
output_time = 500
carb_pr = 5.0 # mm/yr
t = 0
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

V = FunctionSpace(mesh, "CG", 1)

phi_ = Function(V, name="SedOld")
phi = Function(V, name="Sed")
surface = Function(V, name="surface")
limit = Function(V,name="limiter")
thickness = Function(V, name="thickness")
depth = Function(V, name="depth")
nu = Function(V, name="diff")
slf = Function(V, name="sea_level")
light_func = Function(V, name="light_attenuation")
layers = []

v = TestFunction(V)

x = SpatialCoordinate(mesh)
ic = EXAMPLE_INITIAL_COND(x,V)
sl_time = Constant(25 * sin(t / 50000 * 180 / 3.142))
sea_level = Function(V, val=interpolate(sl_time,V), name="sea_level")
# if we don't do this, the sea level ouput has a random function name. 
# Not helpful...
slf.interpolate(sea_level)

land = EXAMPLE_LAND(x, V)
outfile = File("land.pvd")
outfile.write(land)

phi_.assign(ic)
phi.assign(ic)

nu_max = Function(FunctionSpace(mesh, "CG", 1), name='diffusion')
nu_max.assign(1.0)


tiny = 1e-10
surface.interpolate((((land+phi) + land)+abs((land+phi)-land)) / 2)
thickness.interpolate(surface-land)
limit.interpolate((surface - land)/(surface-land+tiny))
depth.interpolate(sea_level - surface)
nu.project(nu_max*((2. / sqrt(2. * math.pi)) * exp(-0.5 * ((depth-5.0)/10.0)**2)))

timestep = 50.0
inflow = 0

F = (inner((phi - phi_)/timestep, v)
     + limit*nu*inner(grad(phi+land), grad(v)))*dx

# We now create an object for output visualisation::

geometry = File("surfaces.pvd")
geometry.write(surface, phi, time=0)
layer_data = File("layer_data.pvd")
layer_data.write(nu, thickness, depth, time=0)
sl = File("Sea_level.pvd")
sl.write(slf, time=0)

end = 50000

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
    functional_expression += inner(phi-pd_gauge_fields[i], phi-pd_gauge_fields[i]) * bump_functions[i] * dx

# Initialise functional to zero
global time_integrated_functional
time_integrated_functional = 0.0


while (t <= run_time):
    # grow carbs
    light_func.interpolate(1.0/(1+exp(-2*depth*25.0))* exp(-1.0*depth/10))
    solve(F == 0, phi)#, bcs=[bc_land, bc_sea])
    #print(phi.dat.data)
    phi += carb_pr * light_func
    #print(phi.dat.data)
    phi_.assign(phi)
    t += timestep
    limit.interpolate((surface - land)/(surface-land+tiny))
    surface.interpolate((((land+phi) + land)+abs((land+phi)-land)) / 2)
    depth.interpolate(sea_level - surface)
    nu.interpolate(nu_max*((2. / sqrt(2. * math.pi)) * exp(-0.5 * ((depth-5.0)/10.0)**2)))
    thickness.interpolate(surface-land)
    sl_time.assign(25 * sin(t / 50000 * 180 / 3.142))
    sea_level.interpolate(sl_time)
    slf.interpolate(sea_level)
    # perform some adjoint magic
    bv = nu.create_block_variable()
    bv.checkpoint = nu._ad_create_checkpoint()
    bv2 = depth.create_block_variable()
    bv2.checkpoint = depth._ad_create_checkpoint()
    bv3 = surface.create_block_variable()
    bv3.checkpoint = surface._ad_create_checkpoint()
    bv4 = limit.create_block_variable()
    bv4.checkpoint = limit._ad_create_checkpoint()
    bv5 = slf.create_block_variable()
    bv5.checkpoint = slf._ad_create_checkpoint()
    bv6 = thickness.create_block_variable()
    bv6.checkpoint = thickness._ad_create_checkpoint()



    if t == 0:
        index = 0
    else:
        index = np.argmax(np.where(times <= t))
    sum_diff = 0
    # update well_from_file
    for i in range(len(coords)):
        pd_gauge_fields[i].assign(wells_from_file[i][index][0])

    time_integrated_functional += assemble(functional_expression)*timestep

    if t%output_time == 0:
        layer_data.write(nu, thickness, depth, time=t)
        geometry.write(surface, phi, time=t)
        sl.write(slf, time=t)
        layers.append([surface, phi, depth, nu, thickness])



time_integrated_functional /= 50000
time_integrated_functional /= timestep

# now we set up the adjoint
c = Control(nu_max)
#djd0 = compute_gradient(time_integrated_functional, c)
#print("Gradient = ",float(djd0))

rf = ReducedFunctional(time_integrated_functional, c)



td_opt = minimise(rf, bounds=[0, 10.0],
                  options={'maxiter': 50})
File('optimal_carb_prod.pvd').write(td_opt)
