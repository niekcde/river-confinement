#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 09:58:07 2024

@author: niekcollotdescury
"""
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
    
    xs, xy = [],[]

    for i in range(samples):

        point = line.interpolate(i  / samples - 1. , normalized=True)
        
        xs.append(point.x)
        xy.append(point.y)
            # access the nearest pixel in the xarray
        tgt_x = xr.DataArray(xs, dims="points")
        tgt_y = xr.DataArray(xy, dims="points")

        dist.append(line.project(point))
    
    profile = xarr.sel(x=tgt_x, y=tgt_y, method="nearest").data
    if len(profile) == 1:
        profile = profile[0]
    
    return dist, profile









