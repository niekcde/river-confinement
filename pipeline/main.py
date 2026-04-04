#!/usr/bin/env python3
# -*- coding: utf-8 -*-
if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

from .build_step2_results import main


if __name__ == "__main__":
    main()
