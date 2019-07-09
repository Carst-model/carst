import os
from typing import Iterable

import firedrake as fd

from .functions import FunctionContainer, carst_funcs

_WANTED_FILES = {
    "land":
    lambda funcs, land: (land, ),
    "layer_data":
    lambda funcs, land: (funcs[carst_funcs.diff_coeff], funcs[
        carst_funcs.thickness], funcs[carst_funcs.depth]),
    "surfaces":
    lambda funcs, land: (funcs[carst_funcs.surface], funcs[carst_funcs.sed]),
    "sea_level":
    lambda funcs, land: (funcs[carst_funcs.sea_level], ),
}


class OutputFilesCollection:
    """Expects a dict for wanted_files
    """

    def __init__(self, output_folder: str,
                 enabled_processes: Iterable[carst_funcs]):
        output_folder = os.path.join(os.getcwd(), output_folder)
        if not os.path.isdir(output_folder):
            raise IOError(f"{output_folder} not a valid path")
        self._out_files = {
            file_name: fd.File(os.path.join(output_folder, file_name + ".pvd"))
            for file_name in _WANTED_FILES
        }

    def output(self,
               funcs: FunctionContainer,
               options,
               names: Iterable[str] = None):
        to_write = self._out_files.keys() if names is None else names
        if not set(to_write).issubset(set(self._out_files.keys())):
            raise AttributeError(
                "Passed names list contains files not managed by this module")

        for file_name in to_write:
            self._out_files[file_name].write(
                *_WANTED_FILES[file_name](funcs, options["land"]),
                time=options["times"]["current_time"],
            )
