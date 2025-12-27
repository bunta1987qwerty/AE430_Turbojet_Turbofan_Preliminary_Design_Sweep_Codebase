# -*- coding: utf-8 -*-
from __future__ import division

import logging
import os
import sys

def setup_logger(out_dir, log_filename='run.log', level=logging.INFO, **kwargs):
    """
    Python 2.7 safe logger creator.

    Accepts both log_filename and legacy keyword names (log_name, etc.) via **kwargs
    so older callers don't crash.
    """
    if not os.path.isdir(out_dir):
        try:
            os.makedirs(out_dir)
        except Exception:
            pass

    # Backward-compat aliases
    if 'log_name' in kwargs and kwargs.get('log_name'):
        log_filename = kwargs.get('log_name')
    if 'logfile' in kwargs and kwargs.get('logfile'):
        log_filename = kwargs.get('logfile')

    logger = logging.getLogger('ae430')
    logger.setLevel(level)
    logger.propagate = False

    # Clear handlers if re-called
    if logger.handlers:
        for h in list(logger.handlers):
            try:
                logger.removeHandler(h)
            except Exception:
                pass

    fmt = logging.Formatter('%(asctime)s | %(levelname)-7s | %(message)s')

    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    log_path = os.path.join(out_dir, log_filename)
    fh = logging.FileHandler(log_path, mode='w')
    fh.setLevel(level)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    logger.info("Logging to %s", log_path)
    return logger
