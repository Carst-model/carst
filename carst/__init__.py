from __future__ import absolute_import
from carst.utility import *
import carst.solver as solver  # NOQA
from carst.callback import DetectorsCallback # NOQA
from carst.functions import FunctionContainer    # NOQA
from carst.functions import carst_funcs as f
import carst.output as output  # NOQA
from carst._version import get_versions
from carst.options import CarstOptions  # NOQA
from carst.processes import (DIFFUSION_EQUATION_GENERIC, INIT_INTERPOLATION_ORDER,
                        advance_carbonates, advance_diffusion)
from carst.solver import CarstModel
import carst.optimisation as optimisation
import carst.exporter as exporter
import os  # NOQA
import datetime  # NOQA
import copy
import math

__version__ = get_versions()['version']
del get_versions

parameters['pyop2_options']['lazy_evaluation'] = False

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
