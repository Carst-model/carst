#!/usr/bin/env python3
import copy
import math

import firedrake as fd
from carst.options import CarstOptions, initialisation_method
from carst.solver import CarstModel

START_TIME = 0
OUTPUT_TIME = 500
TIME_STEP = 50
OUTPUT_FOLDER = "output"


# Example land function generator
def EXAMPLE_LAND(coordinate_space, function_space):
    return fd.project(100 * fd.tanh(0.0005 * (coordinate_space[0] - 6000)),
                      function_space,
                      name="starting_topo")


# Example initial condition
def EXAMPLE_INITIAL_COND(coordinate_space, function_space):
    return fd.project((10000 / (fd.sqrt(2 * math.pi * 250**2)) * fd.exp(-(
        (coordinate_space[0] - 6000)**2) / (2 * 250**2))) +
                      (25000 / (fd.sqrt(2 * math.pi * 1000**2)) * fd.exp(-(
                          (coordinate_space[0] - 4000)**2) / (2 * 1000**2))),
                      function_space)


# Initialise a solver and add land
my_options = CarstOptions(
    initialisation_method.raw_values,
    fd.RectangleMesh(50, 25, 10000, 5000),
    EXAMPLE_LAND,
    "25 * fd.sin(T / 50000 * 180 / math.pi)",
    (
        START_TIME,
        TIME_STEP,
        OUTPUT_TIME,
    ),
    output_folder=OUTPUT_FOLDER,
    diffusion=True,
    carbonates=False,
    diff_coeff=1.0,
)
my_solver_real_scale = CarstModel(my_options)

# DiffuseSolver doesn't play nice with deepcopy apparently
# my_solver_carbonates = copy.deepcopy(my_solver_real_scale)

# Run with a sample initial condition
my_solver_real_scale.set_condition(
    EXAMPLE_INITIAL_COND(my_solver_real_scale.coordinate_space,
                         my_solver_real_scale.function_space))

# Iterate
while my_solver_real_scale.times["current_time"] <= 20000:
    my_solver_real_scale.advance()
