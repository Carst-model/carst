import firedrake as fd
import os
from functions import carst_funcs


_WANTED_FILES = {
    "land": lambda solver: tuple((
        solver.land,
    )),
    "layer_data": lambda solver: tuple((
        solver.funcs[carst_funcs.diff_coeff],
        solver.funcs[carst_funcs.thickness],
        solver.funcs[carst_funcs.depth]
    )),
    "surfaces": lambda solver: tuple((
        solver.funcs[carst_funcs.surface],
        solver.funcs[carst_funcs.sed]
    )),
    "sea_level": lambda solver: tuple((
        solver.funcs[carst_funcs.sea_level],
    )),
}


class OutputFilesCollection:
    """Expects a dict for wanted_files
    """
    def __init__(self, output_folder, enabled_processes):
        if not os.path.isdir(os.getcwd() + "/" + output_folder):
            raise IOError("That folder does not exist!")
        self._out_files = {
            file_name: fd.File(
                f"{output_folder}/{file_name}.pvd"
            ) for file_name in _WANTED_FILES
        }

    def output(self, solver, names):
        if not set(names).issubset(set(self._out_files.keys())):
            raise AttributeError("Passed names list contains files not managed by this module")

        to_write = self._out_files.keys() if names is None else names
        for file_name in to_write:
            self._out_files[file_name].write(
                *_WANTED_FILES[file_name](solver),
                time=solver.current_time,
            )
