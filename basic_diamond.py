#!/usr/bin/env python3
import copy
import math

import firedrake as fd
from carst.options import CarstOptions, initialisation_method
from carst.solver import CarstModel

# Initialise a solver and add land
my_options = CarstOptions(initialisation_method.diamond_default,
                          file="diamond_input.xml")
my_solver_real_scale = CarstModel(my_options)

# Iterate
while my_solver_real_scale.times["current_time"] <= 20000:
    my_solver_real_scale.advance()
