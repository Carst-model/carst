"""
Defines custom callback functions used to compute various metrics at runtime.

"""
from __future__ import absolute_import
from .utility import *
from abc import ABC, abstractproperty, abstractmethod
import h5py
from collections import defaultdict
from .firedrake import *
import numpy as np


class CallbackManager(defaultdict):
    """
    Stores callbacks in different categories and provides methods for
    evaluating them.

    Create callbacks and register them under ``'export'`` mode

    .. code-block:: python

        cb1 = VolumeConservation3DCallback(...)
        cb2 = TracerMassConservationCallback(...)
        cm = CallbackManager()
        cm.add(cb1, 'export')
        cm.add(cb2, 'export')

    Evaluate callbacks, calls :func:`evaluate` method of all callbacks
    registered in the given mode.

    .. code-block:: python

        cm.evaluate('export')

    """
    def __init__(self):
        super(CallbackManager, self).__init__(OrderedDict)

    def add(self, callback, mode):
        """
        Add a callback under the given mode

        :arg callback: a :class:`.DiagnosticCallback` object
        :arg str mode: register callback under this mode
        """
        key = callback.name
        self[mode][key] = callback

    def evaluate(self, mode, index=None):
        """
        Evaluate all callbacks registered under the given mode

        :arg str mode: evaluate all callbacks under this mode
        :kwarg int index: if provided, sets the export index. Default behavior
            is to append to the file or stream.
        """
        for key in sorted(self[mode]):
            self[mode][key].evaluate(index=index)


class DiagnosticHDF5(object):
    """
    A HDF5 file for storing diagnostic time series arrays.
    """
    def __init__(self, filename, varnames, array_dim=1, attrs=None,
                 comm=COMM_WORLD, new_file=True, dtype='d',
                 include_time=True):
        """
        :arg str filename: Full filename of the HDF5 file.
        :arg varnames: List of variable names that the diagnostic callback
            provides
        :kwarg int array_dim: Dimension of the output array. 1 for scalars.
        :kwarg dict attrs: Additional attributes to be saved in the hdf5 file.
        :kwarg comm: MPI communicator
        :kwarg bool new_file: Define whether to create a new hdf5 file or
            append to an existing one (if any)
        """
        self.comm = comm
        self.filename = filename
        self.varnames = varnames
        self.nvars = len(varnames)
        self.array_dim = array_dim
        self.include_time = include_time
        if comm.rank == 0 and new_file:
            # create empty file with correct datasets
            with h5py.File(filename, 'w') as hdf5file:
                if include_time:
                    hdf5file.create_dataset('time', (0, 1),
                                            maxshape=(None, 1), dtype=dtype)
                for var in self.varnames:
                    hdf5file.create_dataset(var, (0, array_dim),
                                            maxshape=(None, array_dim), dtype=dtype)
                if attrs is not None:
                    hdf5file.attrs.update(attrs)

    def _expand_array(self, hdf5file, varname):
        """Expands array varname by 1 entry"""
        arr = hdf5file[varname]
        shape = arr.shape
        arr.resize((shape[0] + 1, shape[1]))

    def _expand(self, hdf5file):
        """Expands data arrays by 1 entry"""
        # TODO is there an easier way for doing this?
        for var in self.varnames:
            self._expand_array(hdf5file, var)
        if self.include_time:
            self._expand_array(hdf5file, 'time')

    def _nentries(self, hdf5file):
        return hdf5file[self.varnames[0]].shape[0]

    def export(self, variables, time=None, index=None):
        """
        Appends a new entry of (time, variables) to the file.

        The HDF5 is updated immediately.

        :arg time: time stamp of entry
        :type time: float
        :arg variables: values of entry
        :type variables: tuple of float
        :kwarg int index: If provided, defines the time index in the file
        """
        if self.comm.rank == 0:
            with h5py.File(self.filename, 'a') as hdf5file:
                if index is not None:
                    nentries = self._nentries(hdf5file)
                    assert index <= nentries, 'time index out of range {:} <= {:} \n  in file {:}'.format(index, nentries, self.filename)
                    expand_required = index == nentries
                    ix = index
                if index is None or expand_required:
                    self._expand(hdf5file)
                    ix = self._nentries(hdf5file) - 1
                if self.include_time:
                    assert time is not None, 'time should be provided as 2nd argument to export()'
                    hdf5file['time'][ix] = time
                for i in range(self.nvars):
                    hdf5file[self.varnames[i]][ix, :] = variables[i]
                hdf5file.close()


class DiagnosticCallback(ABC):
    """
    A base class for all Callback classes
    """

    def __init__(self, solver_obj, array_dim=1, attrs=None,
                 outputdir=None,
                 export_to_hdf5=True,
                 append_to_log=True,
                 include_time=True,
                 hdf5_dtype='d'):
        """
        :arg solver_obj: Thetis solver object
        :kwarg str outputdir: Custom directory where hdf5 files will be stored.
            By default solver's output directory is used.
        :kwarg int array_dim: Dimension of the output array. 1 for scalars.
        :kwarg dict attrs: Additional attributes to be saved in the hdf5 file.
        :kwarg bool export_to_hdf5: If True, diagnostics will be stored in hdf5
            format
        :kwarg bool append_to_log: If True, callback output messages will be
            printed in log
        :kwarg bool include_time: whether to include time in the hdf5 file
        :kwarg hdf5_dtype: Precision to use in hdf5 output: `d` for double
            precision (default), and `f` for single precision
        """
        self.solver_obj = solver_obj
        if outputdir is None:
            self.outputdir = self.solver_obj._options['output_folder']
        else:
            self.outputdir = outputdir
        self.array_dim = array_dim
        self.attrs = attrs
        self.append_to_hdf5 = export_to_hdf5
        self.append_to_log = append_to_log
        self.hdf5_dtype = hdf5_dtype
        self.include_time = include_time
        self._create_new_file = True
        self._hdf5_initialized = False

    def set_write_mode(self, mode):
        """
        Define whether to create a new hdf5 file or append to an existing one

        :arg str mode: Either 'create' (default) or 'append'
        """
        assert mode in ['create', 'append']
        self._create_new_file = mode == 'create'

    def _create_hdf5_file(self):
        """
        Creates an empty hdf5 file with correct datasets.
        """
        if self.append_to_hdf5:
            comm = self.solver_obj.comm
            create_directory(self.outputdir, comm=comm)
            fname = 'diagnostic_{:}.hdf5'.format(self.name.replace(' ', '_'))
            fname = os.path.join(self.outputdir, fname)
            self.hdf_exporter = DiagnosticHDF5(fname, self.variable_names,
                                               array_dim=self.array_dim,
                                               new_file=self._create_new_file,
                                               attrs=self.attrs,
                                               comm=comm, dtype=self.hdf5_dtype,
                                               include_time=self.include_time)
        self._hdf5_initialized = True

    @abstractproperty
    def name(self):
        """The name of the diagnostic"""
        pass

    @abstractproperty
    def variable_names(self):
        """Names of all scalar values"""
        pass

    @abstractmethod
    def __call__(self):
        """
        Evaluate the diagnostic value.

        .. note::
            This method must implement all MPI reduction operations (if any).
        """
        pass

    @abstractmethod
    def message_str(self, *args):
        """
        A string representation.

        :arg args: If provided, these will be the return value from
            :meth:`__call__`.
        """
        return "{} diagnostic".format(self.name)

    def push_to_log(self, time, args):
        """
        Push callback status message to log

        :arg time: time stamp of entry
        :arg args: the return value from :meth:`__call__`.
        """
        print_output(self.message_str(*args))

    def push_to_hdf5(self, time, args, index=None):
        """
        Append values to HDF5 file.

        :arg time: time stamp of entry
        :arg args: the return value from :meth:`__call__`.
        """
        if not self._hdf5_initialized:
            self._create_hdf5_file()
        self.hdf_exporter.export(args, time=time, index=index)

    def evaluate(self, index=None):
        """
        Evaluates callback and pushes values to log and hdf file (if enabled)
        """
        values = self.__call__()
        time = self.solver_obj._times["current_time"]
        if self.append_to_log:
            self.push_to_log(time, values)
        if self.append_to_hdf5:
            self.push_to_hdf5(time, values, index=index)


class DetectorsCallback(DiagnosticCallback):
    """
    Callback that evaluates the specified fields at the specified locations
    """
    def __init__(self, solver_obj,
                 detector_locations,
                 field_names,
                 name,
                 detector_names=None,
                 **kwargs):
        """
        :arg solver_obj: Thetis solver object
        :arg detector_locations: List of x, y locations in which fields are to
            be interpolated.
        :arg field_names: List of fields to be interpolated.
        :arg name: Unique name for this callback and its associated set of
            detectors. This determines the name of the output h5 file
            (prefixed with `diagnostic_`).
        :arg detector_names: List of names for each of the detectors. If not
            provided automatic names of the form `detectorNNN` are created
            where NNN is the (0-padded) detector number.
        :arg kwargs: any additional keyword arguments, see
            :class:`.DiagnosticCallback`.
        """
        # printing all detector output to log is probably not a useful default:
        kwargs.setdefault('append_to_log', False)
        self.field_dims = [solver_obj._funcs[field_name].function_space().value_size
                           for field_name in field_names]
        attrs = {
            # use null-padded ascii strings, dtype='U' not supported in hdf5, see http://docs.h5py.org/en/latest/strings.html
            'field_names': np.array(field_names, dtype='S'),
            'field_dims': self.field_dims,
        }
        super().__init__(solver_obj, array_dim=sum(self.field_dims), attrs=attrs, **kwargs)

        ndetectors = len(detector_locations)
        if detector_names is None:
            fill = len(str(ndetectors))
            self.detector_names = ['detector{:0{fill}d}'.format(i, fill=fill) for i in range(ndetectors)]
        else:
            assert ndetectors == len(detector_names), "Different number of detector locations and names"
            self.detector_names = detector_names
        self._variable_names = self.detector_names
        self.detector_locations = detector_locations
        self.field_names = field_names
        self._name = name

    @property
    def name(self):
        return self._name

    @property
    def variable_names(self):
        return self.detector_names

    def _values_per_field(self, values):
        """
        Given all values evaulated in a detector location, return the values per field"""
        i = 0
        result = []
        for dim in self.field_dims:
            result.append(values[i:i+dim])
            i += dim
        return result

    def message_str(self, *args):
        return '\n'.join(
            'In {}: '.format(name) + ', '.join(
                '{}={}'.format(field_name, field_val) for field_name, field_val in zip(self.field_names, self._values_per_field(values)))
            for name, values in zip(self.detector_names, args))

    def _evaluate_field(self, field_name):
        field = self.solver_obj._funcs[field_name]
        return field(self.detector_locations)

    def __call__(self):
        """
        Evaluate all current fields in all detector locations

        Returns a ndetectors x ndims array, where ndims is the sum of the
        dimension of the fields.
        """
        ndetectors = len(self.detector_locations)
        field_vals = []
        for field_name in self.field_names: 
            field_vals.append(np.reshape(self._evaluate_field(field_name), (ndetectors, -1)))

        return np.hstack(field_vals)


class CompareDetectorsCallback(DetectorsCallback):
    """
    Callback that evaluates the specified fields at the specified locations and compares to previous data
    """
    
    def __init__(self, solver_obj,
                 detector_locations,
                 field_names,
                 functional_expression,
                 name,
                 data,
                 detector_names=None,
                 **kwargs):
        """
        :arg solver_obj: Thetis solver object
        :arg detector_locations: List of x, y locations in which fields are to
            be interpolated.
        :arg field_names: List of fields to be interpolated.
        :arg name: Unique name for this callback and its associated set of
            detectors. This determines the name of the output h5 file
            (prefixed with `diagnostic_`).
        :arg data: filename of previous data (h5 diagnostic file
        :arg detector_names: List of names for each of the detectors. If not
            provided automatic names of the form `detectorNNN` are created
            where NNN is the (0-padded) detector number.
        :arg kwargs: any additional keyword arguments, see
            :class:`.DiagnosticCallback`.
        """

        super().__init__(solver_obj, detector_locations, field_names, name, detector_names, **kwargs)

        # read in previous data
        self.wells_from_file = {}
        self.field_names = field_names
        self.functional = 0
        self.solver_obj = solver_obj
        self.functional_expression = functional_expression
        self.dt = solver_obj._times['time_step']
        f = h5py.File(data, 'r')
        times = np.array(f['time'])
        self.times = times.reshape(times.shape[0])
        self.ndets = len(detector_locations)
        fill = len(str(self.ndets))
        for n in range(self.ndets):
            detector_name = 'detector{:0{fill}d}'.format(n, fill=fill)
            self.wells_from_file[detector_name] = np.array(f[detector_name])
            
    @property
    def rms_error(self):
        """The RMS error between current data and the stored data"""
        # grab data at this timestep
        current_time = self.solver_obj._times['current_time']
        if current_time == 0:
            index = 0
        else:
            index = np.argmax(np.where(self.times <= current_time))
        sum_diff = 0
        fill = len(str(self.ndets))
        current_vals = self.__call__()
        for n in range(self.ndets):
            detector_name = 'detector{:0{fill}d}'.format(n, fill=fill)
            data = self.wells_from_file[detector_name][index][0]
            sum_diff += (data - current_vals[n][0])*(data - current_vals[n][0])

        # perform some adjoint magic
        bv = self.solver_obj._funcs[self.field_names[0]].create_block_variable()
        bv.checkpoint = self.solver_obj._funcs[self.field_names[0]]._ad_create_checkpoint()
        
        rms = sum_diff
        self.functional += assemble(self.functional_expression) * self.dt

        return self.functional

