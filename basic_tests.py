#!/usr/bin/env python3
import copy
import math

import firedrake as fd
from carst.options import CarstOptions
from carst.solver import CarstModel

START_TIME = 0
OUTPUT_TIME = 500
TIME_STEP = 50
OUTPUT_FOLDER = "output"

# Initialise a solver and add land
my_options = CarstOptions(
    fd.RectangleMesh(50, 25, 10000, 5000),
    lambda coord_space, function_space: fd.project(100 * fd.tanh(0.0005 * (
        coord_space[0] - 6000)),
             function_space,
             name="starting_topo"),
    "25 * fd.sin(t / 50000 * 180 / 3.142)",
    (
        START_TIME,
        TIME_STEP,
        OUTPUT_TIME,
    ),
    output_folder=OUTPUT_FOLDER,
    diffusion=True,
    carbonates=False,
    diff_coeff = 1.0,
)
my_solver_real_scale = CarstModel(my_options)

# DiffuseSolver doesn't play nice with deepcopy apparently
# my_solver_carbonates = copy.deepcopy(my_solver_real_scale)

# Run with a sample initial condition
my_solver_real_scale.set_condition(
    fd.project(
        #fd.Constant(5.0),
        (20000 *
         (1 /
          (2 * fd.sqrt(2 * math.pi * 250 * 250))) * fd.exp(-(
              (my_solver_real_scale.coordinate_space[0] - 6000) *
              (my_solver_real_scale.coordinate_space[0] - 6000)) /
                                                           (2 * 250 * 250))) +
        (50000 * (1 / (2 * fd.sqrt(2 * math.pi * 1000 * 1000))) *
         fd.exp(-((my_solver_real_scale.coordinate_space[0] - 4000) *
                  (my_solver_real_scale.coordinate_space[0] - 4000)) /
                (2 * 1000 * 1000))), 
     my_solver_real_scale.function_space))
while my_solver_real_scale.times["current_time"] <= 20000:
    my_solver_real_scale.advance()
