#!/usr/bin/env python3
from collections import UserDict
from enum import Enum
from typing import Callable, Tuple
from xml import etree

import firedrake as fd

from .output import OutputFilesCollection
from .processes import PROCESSOR_NEEDED_FUNCS


class initialisation_method(Enum):
    raw_values = 1
    diamond_default = 2


class CarstOptions(UserDict):
    # Dispatch
    def __init__(
            self,
            ini_type: initialisation_method = initialisation_method.raw_values,
            *args,
            **kw_args):
        if ini_type == initialisation_method.raw_values:
            super().__init__(self.raw_values(*args, **kw_args))
            print("\n".join(f"{repr(key)}: {repr(val)}"
                            for key, val in self.items()))
        elif ini_type == initialisation_method.diamond_default:
            schema_tree = etree.ElementTree(kw_args.get("schema_file"))

    def raw_values(self, base_mesh: fd.mesh.MeshGeometry, land: Callable,
                   sea_level_constant: fd.Constant,
                   times: Tuple[float, float, float], output_folder,
                   **kw_args: dict) -> dict:
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
            if vals["carbonate_production"] is None:
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

        # Initialise out_files
        vals["out_files"] = OutputFilesCollection(output_folder,
                                                  vals["enabled_steps"])
        return vals
