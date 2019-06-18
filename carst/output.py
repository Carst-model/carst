import firedrake as fd


class OutputFilesCollection:
    def __init__(self, solver, wanted_files, output_folder):
        self._write_methods = wanted_files
        self._solver = solver
        self._out_files = {
            file_name: fd.File("{0}/{1}.pvd".format(output_folder, file_name)) for file_name in self._write_methods.keys()
        }

    def output(self):
        for file_name in self._out_files.keys():
            self._out_files[file_name].write(
                *self._write_methods[file_name](self._solver),
                time=self.current_time,
            )

    def output_selective(self, names):
        if not set(names).issubset(set(self._out_files.keys())):
            raise AttributeError("Passed names list contains files not managed by this module")

        for file_name in names:
            self._out_files[file_name].write(
                *self._write_methods[file_name](self._solver),
                time=self.current_time,
            )
