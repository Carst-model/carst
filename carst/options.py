#!/usr/bin/env python3
from collections import UserDict
from typing import Callable, Tuple

import firedrake as fd

from .output import OutputFilesCollection
from .processes import DIFFUSION_EQUATION_GENERIC, PROCESSOR_NEEDED_FUNCS


class CarstOptions(UserDict):
    def __init__(self, base_mesh: fd.mesh.MeshGeometry, land: Callable,
                 sea_level_constant: fd.Constant,
                 times: Tuple[float, float, float], output_folder,
                 **kw_args: dict):
        if not isinstance(base_mesh, fd.mesh.MeshGeometry):
            raise TypeError("base_mesh not of type firedrake.Mesh")
        if not isinstance(sea_level_constant, fd.Constant):
            raise TypeError(
                "sea_level_constant not of type firedrake.Constant")

        vals = dict()

        # Store the passed values
        vals["sea_level_constant"] = sea_level_constant
        vals["times"] = dict(
            zip(("current_time", "time_step", "output_time"), times))
        vals["mesh"] = base_mesh

        # Mark the steps in the process we want
        vals["enabled_steps"] = {
            "diffusion": bool(kw_args.get("diffusion")),
            "carbonates": bool(kw_args.get("carbonates")),
        }
        vals["carbonate_production"] = kw_args.get("carbonate_production")
        if vals["enabled_steps"]["carbonates"]:
            if not vals["carbonate_production"]:
                raise AttributeError(
                    "If carbonate modelling is enabled, a value for the carbonate production rate is required"
                )

        # Generate our workspace from the mesh
        vals["coordinate_space"] = fd.SpatialCoordinate(vals["mesh"])
        vals["function_space"] = fd.FunctionSpace(vals["mesh"], "CG", 1)
        vals["test_function"] = fd.TestFunction(vals["function_space"])

        # Get our land, based in our workspace
        vals["land"] = land(vals["coordinate_space"], vals["function_space"])

        # Initialise the funcs we need
        vals["wanted_funcs"] = list(PROCESSOR_NEEDED_FUNCS["basic"])
        for process in vals["enabled_steps"]:
            if vals["enabled_steps"][process]:
                vals["wanted_funcs"].extend(PROCESSOR_NEEDED_FUNCS[process])

        # Initialise _out_files
        vals["out_files"] = OutputFilesCollection(output_folder,
                                                  vals["enabled_steps"])

        super().__init__(vals)
