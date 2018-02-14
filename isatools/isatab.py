
from __future__ import print_function

import logging

try:
    import pandas as pd
    import numpy as np
    import sys
    print("getting pandas from", pd.__path__, file=sys.stderr)
    from .isatab_meta import *
    from .isatab_full import *
except ImportError as e:
    logging.getLogger(__name__).warn("Missing module")
    logging.getLogger(__name__).warn(str(e))
    logging.getLogger(__name__).warn("Some ISA-tools functionality will not be loaded")
    from .isatab_meta import *
