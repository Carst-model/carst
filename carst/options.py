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


def _process_string_lit(target, replacements):
    result = target
    for original, replacement in replacements:
        result = result.replace(original, replacement)
    return "fd.project(" + result + ", self['function_space'])"


class initialisation_method(Enum):
    """Enum for specifying to carst.options.CarstOptions where you want the initialisation data to come from.
    """
    raw_values = 1
    diamond_default = 2


class CarstOptions(UserDict):
    _STRING_LIT_REPLACEMENTS = (
        ("tanh", "fd.tanh"),
        ("sqrt", "fd.sqrt"),
        ("exp", "fd.exp"),
        ("sin", "fd.sin"),
        ("X", "self['coordinate_space'][0]"),
        ("F", "self['function_space']"),
        ("pi", "math.pi"),
    )

    # Dispatch
    def __init__(self, ini_type: initialisation_method, *args, **kw_args):
        """Initialise the carst.options.CarstOptions module.

        If ini_type is carst.options.initialisation_method.raw_values, args should contain:

        * The base mesh on which the model will operate, of type firedrake.mesh.MeshGeometry.

        :param carst.options.initialisation_method ini_type: The method to use for gathering data.
        :param Iterable args: Only used if if ini_type is carst.options.initialisation_method.raw_values.
        :returns: The processed carst.options.CarstModel instance.
        """
        super().__init__()
        self["type"] = ini_type
        if ini_type == initialisation_method.raw_values:
            self._raw_values(*args, **kw_args)
        elif ini_type == initialisation_method.diamond_default:
            self._diamond_default(kw_args["file"])

    def _diamond_default(self, file_name="diamond_input.xml"):
        tree_root = ElementTree.ElementTree(file=file_name).getroot()
        del tree_root[0]

        sea_level_lit = tree_root[5][0].text
        for original, replacement in CarstOptions._STRING_LIT_REPLACEMENTS:
            sea_level_lit = sea_level_lit.replace(original, replacement)
        self["sea_level"] = sea_level_lit.replace(
            "T", "self._times['current_time']")

        make_int = lambda num: int(num)
        mesh_args = list(
            zip(map(make_int, tree_root[0][0][0].text.split(" ")),
                map(make_int, tree_root[0][1][0].text.split(" "))))
        mesh_args = mesh_args[0] + mesh_args[1]
        self["mesh"] = fd.RectangleMesh(*mesh_args)
        self["coordinate_space"] = fd.SpatialCoordinate(self["mesh"])
        self["function_space"] = fd.FunctionSpace(self["mesh"], "CG", 1)
        self["test_function"] = fd.TestFunction(self["function_space"])

        self["times"] = dict(
            zip(("current_time", "time_step", "output_time", "end_time"),
                [int(time[0].text) for time in tree_root[1]]))

        self["enabled_steps"] = {
            step.tag.lower(): step.text == "True"
            for step in tree_root[6]
        }
        print(self["enabled_steps"])
        # Locate diff_coeff in optional parts of the tree
        if self["enabled_steps"]["diffusion"]:
            diff_coeff_not_present = AttributeError(
                "Diffusion enabled but no coefficient provided!")
            try:
                if tree_root[7].tag == "Diff-Coefficient":
                    self["diff_coeff"] = float(tree_root[7][0].text)
                elif tree_root[8].tag == "Diff-Coefficient":
                    self["diff_coeff"] = float(tree_root[8][0].text)
                else:
                    raise diff_coeff_not_present
            except IndexError:
                raise diff_coeff_not_present
        # Initialise the funcs we need
        self["wanted_funcs"] = list(PROCESSOR_NEEDED_FUNCS["basic"])
        for process in self["enabled_steps"]:
            if self["enabled_steps"][process]:
                self["wanted_funcs"].extend(PROCESSOR_NEEDED_FUNCS[process])

        self["out_files"] = OutputFilesCollection(tree_root[2][0].text,
                                                  self["enabled_steps"])

        # Evaluate the literals for the initial_condition and land
        self["land"] = eval(
            _process_string_lit(tree_root[3][0].text,
                                CarstOptions._STRING_LIT_REPLACEMENTS))
        self["initial_condition"] = eval(
            _process_string_lit(tree_root[4][0].text,
                                CarstOptions._STRING_LIT_REPLACEMENTS))

    def _raw_values(self, base_mesh, land: Callable, sea_level,
                    times: Tuple[float, float, float], output_folder,
                    **kw_args: dict) -> dict:
        if not isinstance(base_mesh, fd.mesh.MeshGeometry):
            raise TypeError("base_mesh not of type firedrake.Mesh")
        if not isinstance(sea_level, str):
            raise TypeError("sea_level not of type str")

        # Store the passed values
        self["sea_level"] = sea_level.replace("T",
                                              "self._times['current_time']")
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

        if kw_args.get("initial_condition") is not None:
            self["initial_condition"] = kw_args["initial_condition"]
