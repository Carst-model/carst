from __future__ import absolute_import
from carst.utility import *
from carst.log import *
# Note: NOQA is No Quality Assurance. The linter ignores warnings from those lines
import carst.timeintegrator as timeintegrator  # NOQA
import carst.solver2d as solver2d  # NOQA
from carst.callback import DiagnosticCallback  # NOQA
import carst.limiter as limiter      # NOQA
from carst._version import get_versions
__version__ = get_versions()['version']
del get_versions

parameters['pyop2_options']['lazy_evaluation'] = False

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
