#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 10:38:35 2024

@author: niekcollotdescury
"""

import os
import glob

def remove_file(path):
    if os.path.exists(path):
        os.remove(path)
    else:
        print(f"The file does not exist: {path}")

def remove_dems(df, directory):
    ids = df.reach_id.values
    for id in ids:
        files = glob.glob(f'{directory}input_created/dem/{id}_*')
        for f in files:
            remove_file(f)