# Import firedrake stuff here so that adjoint mode works
from __future__ import absolute_import

import carst_config

from firedrake import *  # NOQA

if carst_config.adjoint:
    from firedrake_adjoint import *  # NOQA

