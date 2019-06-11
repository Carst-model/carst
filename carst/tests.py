#!/usr/bin/env python3
import math
import firedrake as fd
from carst_it1 import Diffuse_Solver

# Initialise a solver and add land
my_solver = Diffuse_Solver(fd.RectangleMesh(50, 25, 10000, 5000), "output")
my_solver.add_land(fd.project(
    100 * fd.tanh(0.0005 * (my_solver.coordinate_space[0] - 6000)),
    my_solver.function_space,
    name="starting_topo",
))

# Run with a sample initial condition
my_solver.diffuse_real_scale_test(
    fd.project(
        (
            20000
            * (1 / (2 * fd.sqrt(2*math.pi*250*250)))
            * fd.exp(
                -((my_solver.coordinate_space[0]-6000) * (my_solver.coordinate_space[0]-6000))
                / (2 * 250 * 250)
            )
        )
        + (
            50000
            * (1 / (2 * fd.sqrt(2 * math.pi * 1000 * 1000)))
            * fd.exp(
                -((my_solver.coordinate_space[0]-4000) * (my_solver.coordinate_space[0]-4000))
                / (2 * 1000 * 1000)
            )
        ), my_solver.function_space),
    0,
    20000,
    50,
    500,
)
