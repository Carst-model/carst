#!/usr/bin/env python3
import math
import copy
import firedrake as fd
from carst import DiffuseSolver

# Initialise a solver and add land
my_solver_real_scale = DiffuseSolver(
    fd.RectangleMesh(50, 25, 10000, 5000),
    "output",
)
my_solver_real_scale.add_land(
    fd.project(
        100 * fd.tanh(0.0005 * (my_solver_real_scale.coordinate_space[0] - 6000)),
        my_solver_real_scale.function_space,
        name="starting_topo",
    ),
)

# DiffuseSolver doesn't play nice with deepcopy apparently
# my_solver_carbonates = copy.deepcopy(my_solver_real_scale)

# Run with a sample initial condition
my_solver_real_scale.diffuse_real_scale_test(
    fd.project(
        (
            20000
            * (1 / (2 * fd.sqrt(2*math.pi*250*250)))
            * fd.exp(
                -((my_solver_real_scale.coordinate_space[0]-6000) * (my_solver_real_scale.coordinate_space[0]-6000))
                / (2 * 250 * 250)
            )
        )
        + (
            50000
            * (1 / (2 * fd.sqrt(2 * math.pi * 1000 * 1000)))
            * fd.exp(
                -((my_solver_real_scale.coordinate_space[0]-4000) * (my_solver_real_scale.coordinate_space[0]-4000))
                / (2 * 1000 * 1000)
            )
        ), my_solver_real_scale.function_space),
    0,
    20000,
    50,
    500,
)

# Run with a sample initial condition
# my_solver_carbonates.diffuse_real_scale_test(
    # fd.project(
        # (
            # 20000
            # * (1 / (2 * fd.sqrt(2*math.pi*250*250)))
            # * fd.exp(
                # -((my_solver_carbonates.coordinate_space[0]-6000) * (my_solver_carbonates.coordinate_space[0]-6000))
                # / (2 * 250 * 250)
            # )
        # )
        # + (
            # 50000
            # * (1 / (2 * fd.sqrt(2 * math.pi * 1000 * 1000)))
            # * fd.exp(
                # -((my_solver_carbonates.coordinate_space[0]-4000) * (my_solver_carbonates.coordinate_space[0]-4000))
                # / (2 * 1000 * 1000)
            # )
        # ), my_solver_carbonates.function_space),
    # 0,
    # 20000,
    # 50,
    # 500,
# )
