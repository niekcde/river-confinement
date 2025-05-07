#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 10:30:53 2024

@author: niekcollotdescury
"""

import scipy as sc
import shapely

from shapely.geometry import LineString


def apply_smoothing(line, w):
    if w > len(line.coords):
        w = len(line.coords)
    if w == 0:
        w = 2
    w = int(w)

    XSG = sc.signal.savgol_filter(line.xy[0], w, 1)
    YSG = sc.signal.savgol_filter(line.xy[1], w, 1)
    SG  = LineString(list(zip(XSG, YSG)))
    
    return SG
import matplotlib.pyplot as plt
def SG_smoothing(line, w, width):
    """Function that used a savtisky Golay filter from scipy to smooth a linestring.
    A hausdorf distance larger than 0.5*width creates new iteration with smaller smoothing window
    input:
        line: shapely LineString
        w: smoothing window
        width: check value for hausdorf distance
    return: Smoothed LineString
    """
    line = line.segmentize(1)
    line       = line.simplify(0.1, preserve_topology=True)

    lineSmooth = apply_smoothing(line, w)

    hdist = shapely.hausdorff_distance(line, lineSmooth)
    while hdist > (width*0.5):
        w *= 0.95        
        lineSmooth = apply_smoothing(line, w)
        
        hdist = shapely.hausdorff_distance(line, lineSmooth)
        
    return lineSmooth



# def SG_smoothing(line, wInn, full = False, lowerBoundary = 0.02,upperBoundary = 0.06):
    """Function that used a savtisky Golay filter from scipy to smooth a linestring
    input:
        line: shapely LineString
        w: smoothing window
        full: if full is True the line will be segmentized to 1 meter segments if False
        lowerboundary: minimum relative size smoothing window (default 2%)
        upperboundary: maximum relative size smoothing window (default 6%)
    return: Smoothed LineString
    """
    #########################
    # full segmentation or use reduced linestring for efficiency
    #########################
    w = wInn

    # if full == False:
    #     line = line.simplify(0.5, preserve_topology=True)

    #     segment_dis = 25
    #     max_segments = line.segmentize(1)
    #     line  = line.segmentize(segment_dis)
    #     ratio = int(len(max_segments.coords) / len(line.coords))

    #     w = int(w / ratio)
    #     if w == 1:
    #         w = 2
            


    # else:
    #     line  = line.segmentize(2)

    line = line.segmentize(1)

    #########################
    # limit the smoothing window to upper and lower boundary
    #########################
    # coords = len(line.coords)
    # count, countif1, countif2 = 0,0,0

    # while ( ( (w / coords) < lowerBoundary) | ((w / coords) > upperBoundary)):
    #     if (w / coords) < lowerBoundary:
    #         w = w*1.3
    #         countif1 += 1
    #     if (w / coords) > upperBoundary:
    #         w = w *0.9
    #         countif2 += 1

    #     if (countif1 > 5) & (countif2 > 5):
    #         w = (lowerBoundary + upperBoundary) /2
    #         break
    #     count += 1
    #     if count == 1e4:
    #         break
    #########################
    # Check: if window size to large convert to max coords
    #########################
    if w > len(line.coords):
        w = len(line.coords)
    if w == 0:
        w = 2
    
    #########################
    # Apply savitzky golay filter
    #########################
    w = int(w)
    XSG = sc.signal.savgol_filter(line.xy[0], w, 1)
    YSG = sc.signal.savgol_filter(line.xy[1], w, 1)
    SG  = LineString(list(zip(XSG, YSG)))
    
    line = SG.simplify(0.01, preserve_topology=True)
    # line = SG
    return line
