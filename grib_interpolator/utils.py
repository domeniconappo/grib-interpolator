"""
This software comes as Open Source and licensed via AGPL v3.
It was developed under the initiative Copernicus, EFAS operational center @ECMWF (Reading, UK).
"""

import numpy as np
from datetime import datetime

int_fill_value = -999999

ENDC = '\033[0m'
BOLD = '\033[1m'
GREEN = '\033[92m' + BOLD
YELLOW = '\033[93m' + BOLD
WARN = '\033[94m' + BOLD
FAIL = '\033[91m' + BOLD
DEFAULT = BOLD


def is_container(a):
    if isinstance(a, (list, tuple, dict)):
        return True
    return False


def is_callable(v):
    return hasattr(v, '__call__')


def mask_it(v, mv, shape=None):
    if shape is not None:
        result = np.ma.masked_array(data=v, fill_value=mv, copy=False)
    else:
        result = np.ma.masked_values(v, mv, copy=False)
    return result


def empty(shape, fill_value=np.NaN, dtype=float):
    idxs = np.empty(shape, dtype=dtype)
    idxs.fill(fill_value)
    return idxs


def now_string(fmt='%Y-%m-%d %H:%M'):
    return datetime.strftime(datetime.now(), fmt)


def progress_step_and_backchar(num_cells):
    progress_step = num_cells / 250
    back_char = '\r'
    # if not stdout.isatty():
    #     # out is being redirected
    #     back_char = '\n'
    #     progress_step *= 10
    return back_char, progress_step

