#!/usr/bin/env python3
import math
from collections import UserDict
from enum import Enum
from itertools import chain
from typing import Callable, Tuple
from xml.etree import ElementTree

import firedrake as fd

from .output import OutputFilesCollection
from .processes import PROCESSOR_NEEDED_FUNCS


class initialisation_method(Enum):
    raw_values = 1
    diamond_default = 2


class CarstOptions(UserDict):
    _STRING_LIT_REPLACEMENTS = (
        ("tanh", "fd.tanh"),
        ("sqrt", "fd.sqrt"),
        ("exp", "fd.exp"),
        ("x", "self['coordinate_space'][0]"),
        ("pi", "math.pi"),
    )

    # Dispatch
    def __init__(
            self,
            ini_type: initialisation_method = initialisation_method.raw_values,
            *args,
            **kw_args):
        super().__init__()
        self["type"] = ini_type
        if ini_type == initialisation_method.raw_values:
            self._raw_values(*args, **kw_args)
        elif ini_type == initialisation_method.diamond_default:
            self._diamond_default(kw_args["file"])

    def _diamond_default(self, file_name):
        tree_root = ElementTree.ElementTree(file="diamond_input.xml").getroot()
        del tree_root[0]

        self["sea_level"] = tree_root[5][0].text
        self["mesh"] = fd.RectangleMesh(*chain(
            zip([int(num) for num in tree_root[0][0][0].text.split(" ")],
                [int(num) for num in tree_root[0][1][0].text.split(" ")])))
        self["function_space"] = fd.FunctionSpace(self["mesh"], "CG", 1)
        self["test_function"] = fd.TestFunction(self["function_space"])

        self["times"] = dict(
            zip(("current_time", "time_step", "output_time"),
                [int(time[0].text) for time in tree_root[1]]))

        # Enabled-Steps must always be the last element of the xml tree
        self["enabled_steps"] = {
            step.tag.lower(): bool(step.text)
            for step in tree_root[-1]
        }

        self["out_files"] = OutputFilesCollection(tree_root[2][0].text,
                                                  self["enabled_steps"])

        land_func = tree_root[6][0].text
        for original, replacement in CarstOptions._STRING_LIT_REPLACEMENTS:
            land_func = land_func.replace(original, replacement)
        land_func = "fd.project(" + land_func + ")"
        self["land"] = eval(land_func)

    def _raw_values(self, base_mesh: fd.mesh.MeshGeometry, land: Callable,
                    sea_level: fd.Constant, times: Tuple[float, float, float],
                    output_folder, **kw_args: dict) -> dict:
        if not isinstance(base_mesh, fd.mesh.MeshGeometry):
            raise TypeError("base_mesh not of type firedrake.Mesh")
        if not isinstance(sea_level, str):
            raise TypeError("sea_level not of type str")

        self = dict()

        # Store the passed values
        self["sea_level"] = sea_level
        self["times"] = dict(
            zip(("current_time", "time_step", "output_time"), times))
        self["mesh"] = base_mesh

        # Mark the steps in the process we want
        self["enabled_steps"] = {
            "diffusion": bool(kw_args.get("diffusion")),
            "carbonates": bool(kw_args.get("carbonates")),
        }
        self["carbonate_production"] = kw_args.get("carbonate_production")
        if self["enabled_steps"]["carbonates"]:
            if self["carbonate_production"] is None:
                raise AttributeError(
                    "If carbonate modelling is enabled, a value for the carbonate production rate is required"
                )

        # Generate our workspace from the mesh
        self["coordinate_space"] = fd.SpatialCoordinate(self["mesh"])
        self["function_space"] = fd.FunctionSpace(self["mesh"], "CG", 1)
        self["test_function"] = fd.TestFunction(self["function_space"])

        # Get our land, based in our workspace
        self["land"] = land(self["coordinate_space"], self["function_space"])

        # Initialise the funcs we need
        self["wanted_funcs"] = list(PROCESSOR_NEEDED_FUNCS["basic"])
        for process in self["enabled_steps"]:
            if self["enabled_steps"][process]:
                self["wanted_funcs"].extend(PROCESSOR_NEEDED_FUNCS[process])

        # Initialise out_files
        self["out_files"] = OutputFilesCollection(output_folder,
                                                  self["enabled_steps"])
        self["diff_coeff"] = float(kw_args.get('diff_coeff'))
