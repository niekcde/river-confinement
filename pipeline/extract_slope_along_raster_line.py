#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 09:58:07 2024

@author: niekcollotdescury
"""
import numpy as np
import shapely
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
    sample_positions = np.arange(samples, dtype=float) / samples - 1.0
    points = shapely.line_interpolate_point(line, sample_positions, normalized=True)

    xs = shapely.get_x(points)
    ys = shapely.get_y(points)
    tgt_x = xr.DataArray(xs, dims="points")
    tgt_y = xr.DataArray(ys, dims="points")
    # Preserve projection semantics for self-intersections and the list return type.
    dist = shapely.line_locate_point(line, points).tolist()

    profile = xarr.sel(x=tgt_x, y=tgt_y, method="nearest").data
    if len(profile) == 1:
        profile = profile[0]
    
    return dist, profile







