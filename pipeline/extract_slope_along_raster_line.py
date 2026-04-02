#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 09:58:07 2024

@author: niekcollotdescury
"""
import numpy as np
import xarray as xr


def extract_slope_along_raster_line(xarr, line, samples = 400):
    ''''Function for extracting a slope along a line segment from a digital elevation Slope raster
    Input:
    - xarr: DEM raster
    - line: Line segment for slope
    - samples: number of sample points for the slope line (200)\n
    return:
    - Distance values for the slope line
    - Elevation values for the slope Line
        '''
    profile = []
    dist    = []

    sample_positions = np.arange(samples, dtype=float) / samples - 1.0
    points = [line.interpolate(position, normalized=True) for position in sample_positions]

    xs = [point.x for point in points]
    ys = [point.y for point in points]
    tgt_x = xr.DataArray(xs, dims="points")
    tgt_y = xr.DataArray(ys, dims="points")
    dist = [line.project(point) for point in points]

    profile = xarr.sel(x=tgt_x, y=tgt_y, method="nearest").data
    if len(profile) == 1:
        profile = profile[0]
    
    return dist, profile








