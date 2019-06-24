#!/usr/bin/env python3
import firedrake as fd
from functions import FunctionContainer
from output import OutputFilesCollection
from processors import PROCESSOR_NEEDED_FUNCS


class CarstOptions:
    def __init__(self, base_mesh, land, sea_level_constant, times, **kw_args):
        if not isinstance(base_mesh, fd.mesh.MeshGeometry):
            raise TypeError("base_mesh not of type firedrake.Mesh")
        if not isinstance(sea_level_constant, fd.Constant):
            raise TypeError("sea_level_constant not of type firedrake.Constant")

        # Store the passed values
        self._sea_level_constant = sea_level_constant
        self.current_time, self._time_step, self._output_time = times
        self.mesh = base_mesh

        # Mark the steps in the process we want
        self.enabled_steps = {
            "diffusion": bool(kw_args.get("diffusion")),
            "carbonates": bool(kw_args.get("carbonates")),
        }
        if self.enabled_steps["carbonates"]:
            self._carbonate_production = kw_args.get("carbonate_production")
            if self._carbonate_production is None:
                raise AttributeError("If carbonate modelling is enabled, a value for the carbonate production rate is required")
        else:
            self._carbonate_production = None

        # Generate our workspace from the mesh
        self.coordinate_space = fd.SpatialCoordinate(self.mesh)
        self.function_space = fd.FunctionSpace(self.mesh, "CG", 1)
        self.test_function = fd.TestFunction(self.function_space)

        # Get our land, based in our workspace
        self._land = land(self.coordinate_space, self.function_space)

        # Initialise the funcs we need
        self._funcs = FunctionContainer(
            self,
            {
                func_name for processor in PROCESSOR_NEEDED_FUNCS.values() for func_name in processor
            },
        )

        # Initialise _out_files
        if kw_args.get("output_folder") is not None:
            self._out_files = OutputFilesCollection(kw_args.get("output_folder"), self.enabled_steps)

    @property
    def useful_data(self):
        return {
            "workspace": tuple((   # Variables for the function space + constants
                self.mesh,
                self.coordinate_space,
                self.function_space,
                self.test_function,
                self._sea_level_constant,
                self._land,
            )),
            "times": tuple((   # The timing-related vars
                self.current_time,
                self._output_time,
                self._time_step,
            )),
            "funcs": self._funcs,
            "enabled_steps": self.enabled_steps,
            "output_files": self._out_files,
            "optional": tuple((
                self._carbonate_production,
            )),
        }
