#!/usr/bin/env python3
import copy
import math
import random

import firedrake as fd
from carst.options import CarstOptions, initialisation_method
from carst.solver import CarstModel

START_TIME = 0
OUTPUT_TIME = 500
TIME_STEP = 50
END_TIME = 10000
OUTPUT_FOLDER = "output"


# Example land function generator
def EXAMPLE_LAND(coordinate_space, function_space):
    return fd.project(100 * fd.tanh(0.0005 * (coordinate_space[0] - 6000) + random.randint(0,10)/10. ),
                      function_space,
                      name="starting_topo")


# Example initial condition
def EXAMPLE_INITIAL_COND(coordinate_space, function_space):
    return fd.project(
        (20000 * (1 / (2 * fd.sqrt(2 * math.pi * 250 * 250))) * fd.exp(-(
            (coordinate_space[0] - 6000) *
            (coordinate_space[0] - 6000)) / (2 * 250 * 250))) +
        (50000 * (1 / (2 * fd.sqrt(2 * math.pi * 1000 * 1000))) * fd.exp(-(
            (coordinate_space[0] - 4000) *
            (coordinate_space[0] - 4000)) / (2 * 1000 * 1000))),
        function_space)


# Initialise a solver and add land
my_options = CarstOptions(
    initialisation_method.raw_values,
    fd.RectangleMesh(200, 200, 10000, 10000),
    EXAMPLE_LAND,
    "25 * fd.sin(t / 50000 * 180 / 3.142)",
    (
        START_TIME,
        TIME_STEP,
        OUTPUT_TIME,
        END_TIME,
    ),
    output_folder=OUTPUT_FOLDER,
    diffusion=True,
    carbonates=True,
    diff_coeff=5.0,
    carbonate_production = 3.0,
)
my_solver_real_scale = CarstModel(my_options)

# DiffuseSolver doesn't play nice with deepcopy apparently
# my_solver_carbonates = copy.deepcopy(my_solver_real_scale)

# Run with a sample initial condition
my_solver_real_scale.set_condition(
    EXAMPLE_INITIAL_COND(my_solver_real_scale.coordinate_space,
                         my_solver_real_scale.function_space))

my_solver_real_scale.iterate()
