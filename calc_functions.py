#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 13:47:48 2024

@author: niekcollotdescury
"""

import numpy as np
import math as m

def z_score(x):
    return (x - np.mean(x)) / np.std(x)

def minmax(x):
    return (x - np.min(x)) / (np.max(x) - np.min(x))

def slope_with_intercept(x:'np.array',y:'np.array', intercept:'float'):
    """
    Calculate the slope of a line with a given intercept.\n
    input:
    - x: x-axis
    - y: y-axis
    - intercept: y axis intercept\n
    return: slope value units of meter over meter

    ADD CHECK FOR CORRECT types
    """
    # Adjust y values by subtracting the intercept
    y_adjusted = y - intercept

    # Calculate the best-fit slope
    slope = np.dot(x, y_adjusted) / np.dot(x, x)
    return(slope)

def get_entrenchment_slope_intersect(distance, y0,y1, centerInt, xintercept, side):
    if side == 'positive':
        # xind = np.where(distance <= xintercept)[0][-1]
        # x0 = distance[centerInt]
        # x1 = distance[xind + 1]

        Xdist = xintercept - distance[centerInt]
        Ydiff = y1 - y0

       

    elif side == 'negative':
        Xdist = abs(xintercept) - abs(distance[centerInt])
        Ydiff = y1 - y0


    slope = Ydiff / Xdist


    return np.rad2deg(np.arctan(slope)) #units of degree

def get_entrenchment_slope_no_intersect(distance, profile, centerDistance, xintercept, yintercept, side):
    distance = distance.copy()
    profile  = profile.copy()
    
    if side == 'positive':
        X = distance[centerDistance ::]
        Y = profile[centerDistance ::]
        Yintercept = Y[0]

    elif side == 'negative':
        X = distance[0:centerDistance+1]
        Y = profile[0:centerDistance+1]

        X = X *-1
        Yintercept = Y[-1]


    slope = slope_with_intercept(X, Y, Yintercept)
    return np.rad2deg(np.arctan(slope))

def linear_fit(x, slope, intercept):
    X = np.array(x)

    return slope * X + intercept

def slope_curvature(xInn, yInn, yTop, xCentInt, xMax):

    if xMax < 0:
        xind = np.where(xInn >= xMax)[0][0]
        if xind > 0:
            xind -= 1
        x    = xInn[xind : xCentInt +1]
        y    = yInn[xind : xCentInt +1]
        
        x = np.flip(x*-1, 0)
        y = np.flip(y, 0)
        xMax = abs(xMax)
    else:
        xind = np.where(xInn <= xMax)[0][-1]
        x = xInn[xCentInt:xind+1]
        y = yInn[xCentInt:xind+1]


    
    y[-1] = yTop
    x[-1] = xMax


    yNorm = y - y.min() # change y values to values with min of zero
    area_under = np.trapz(y = yNorm, x = x)
    area_total = (x.max() - x.min()) * (yNorm.max())


    res = (area_total - area_under) / (area_total / 2)
    return res

def length_averaged_section(vals, sections,line, points =True):

    """
    function that calculates a weighted average
    input:
        vals: values to be averaged
        sections: sections lengths or section points
        line: centerline
        points: determines if sections is supplied as points or as lengths (default = True, points)
    return:
        weighted average of vals
    """
    if points == True:
        sec_lens = []
        for i in range(len(sections) -1):
            dist = abs(line.project(sections[i+1]) - line.project(sections[i]))
            sec_lens.append(dist)
        sec_lens = np.array(sec_lens)
    else:
        sec_lens = sections

    weighted_average = np.average(vals, weights=sec_lens)
    return weighted_average

def calculate_point_rightsided_triangle(x_a, y_a, x_b, y_b, length_AC):
    """
    Calculate the coordinates of point C in a right triangle.\n
    
    Args:
    - x_a, y_a: Coordinates of point A (right angle).
    - x_b, y_b: Coordinates of point B.
    - length_AC: Length of side AC.\n

    Returns:
    - (x_c1, y_c1), (x_c2, y_c2): Two possible coordinates for point C.
    """
    # Calculate the vector AB
    dx = x_b - x_a
    dy = y_b - y_a
    
    # Normalize vector AB
    distance_AB = m.sqrt(dx**2 + dy**2)
    if distance_AB == 0:
        raise ValueError("Points A and B cannot be the same.")
    
    unit_dx = dx / distance_AB
    unit_dy = dy / distance_AB
    
    # Calculate the perpendicular vector to AB
    perp_dx = -unit_dy
    perp_dy = unit_dx

    # Calculate the two possible positions for point C
    x_c1 = x_a + length_AC * perp_dx
    y_c1 = y_a + length_AC * perp_dy
    
    x_c2 = x_a - length_AC * perp_dx
    y_c2 = y_a - length_AC * perp_dy
    
    return (x_c1, y_c1), (x_c2, y_c2)

def extend_apex(apex , apexOrig, inf, width, amplitude, distance):
    """
    Extend Apex in orthogonal directions. Limit distance to 30 kilometer.\n
    input: 
    - apex: apex point
    - apexOrig: Origin point on the inflection line
    - inf: inflection point
    - width: river width value
    - distance: distance factor\n
    output:
    - D: extension outside bend
    - E: extension innside bend
    """
    extendDistance = width*distance
    if extendDistance > 30000:
        extendDistance = 30000

    apexX, apexY         = apex.x    , apex.y
    apexOrigX, apexOrigY = apexOrig.x, apexOrig.y
    infX, infY           = inf.xy[0][0]     , inf.xy[0][1]

    if (amplitude == 0):
        # Calculate the vector AB
        dx = infX - apexX
        dy = infY - apexY
        
        # Normalize vector AB
        distance_AB = m.sqrt(dx**2 + dy**2)
        if distance_AB == 0:
            raise ValueError("Points A and B cannot be the same.")
        
        unit_dx = dx / distance_AB
        unit_dy = dy / distance_AB
        
        # Calculate the perpendicular vector to AB
        perp_dx = -unit_dy
        perp_dy = unit_dx
        # Calculate the two possible positions for point C
        D = apexX + extendDistance * perp_dx, apexY + extendDistance * perp_dy
        E = apexX - extendDistance * perp_dx, apexY - extendDistance * perp_dy
    else:
        dx = apexOrigX - apexX
        dy = apexOrigY - apexY
        d = m.sqrt(dx**2 + dy**2)
        
        D = apexX - (extendDistance * (dx/d)), apexY - (extendDistance * (dy/d))
        E = apexX + (extendDistance * (dx/d)), apexY + (extendDistance * (dy/d))

    return D, E

def x_y_intercept(x, y , threshold, w, adjust_y = False, adjust_height = 0):
    """
    Calculate the intecept between x and y, for y intercept value.\n
    input:
    - x: list of values
    - y: list of values
    - threshold: single value corresponding to value in y axis\n
    output:
    - x intercept or nan for no intercept
    The list of x and y should be of equal length
    """
    x, y = np.array(x), np.array(y)

    if adjust_y == True:
        y[x <= w] = adjust_height

    if np.max(y) > threshold:
        xp = np.argwhere(y > threshold)[0]
        ind = xp[0]
        if ind == 0:
            x_t = w
        x_t = x[ind-1] + (((threshold - y[ind-1]) * (x[ind] - x[ind-1])) / (y[ind] - y[ind-1]))
        return x_t
    else:
        return np.nan

def confinement_slope(intercept, profile, distance,centerPointHeight, confHeight, widthT, normalize = True):
    heightDiff = confHeight - centerPointHeight
    maxSlope = heightDiff / (widthT / 2)

    if np.isnan(intercept):
        slope = slope_with_intercept(distance, profile, centerPointHeight)
    else:
        
        
        if intercept < (widthT / 2):
            slope = maxSlope
        else:
            slope = heightDiff / intercept

    if normalize == True:
        normSlope = slope / maxSlope
    else:
        normSlope = slope
    return normSlope

def confinement_ratio(intercept, width):
    if np.isnan(intercept):
        ratio = np.nan
    else:
        ratio = intercept / (width)
    return ratio

def get_low_center_point(po, pi, cdo, cdi, width):
    if (isinstance(po, int)) & (isinstance(pi, int)):
        return np.nan
    elif (isinstance(po, list)) & (isinstance(pi, int)):
        po  = np.array(po)
        cdo = np.array(cdo)
        poRiv = po[cdo <= (width / 2)]
        return np.min(poRiv)
    elif (isinstance(po, int)) & (isinstance(pi, list)):
        pi  = np.asarray(pi)
        cdi = np.asarray(cdi)
        piRiv = pi[cdi <= (width / 2)]
        return np.min(piRiv)
    else:

        po, pi   = np.array(po), np.asarray(pi)
        cdo, cdi = np.array(cdo), np.asarray(cdi)

        poRiv = po[cdo <= (width / 2)]
        piRiv = pi[cdi <= (width / 2)]

        return np.min([poRiv.min(), piRiv.min()])

def confinement_values(po, pi, cdo, cdi,widthW, widthT, factor):
    """
    Calculate confinement slope and ratio of intercept with river width.\n
    input:
    - po/pi: elevation profile for outer and inner bend
    - cdo/cdi: distance from river centerline for outer and inner bend
    - width: river width
    - factor: factor that multiplied with the river width determines the confinement height\n
    output:
    - po/pi Intercept: intercept distance for outer and inner bend
    - centerPointHeight: height of the centerpoint of the river determmined by the lowest point within river bounds
    - confinementHeight: valley edge height based on factor multiplied with river width above river centerpoint
    - slope Out/Inn: normalized outer and inner slope
    - ER Out/Inn: Ratio of valley edge intercept and river width for outer and inner bend
    """

    centerPointHeight = get_low_center_point(po, pi, cdo, cdi, widthT)
    confinementHeight = (widthW * factor) + centerPointHeight

    if (isinstance(po, int)):
        poIntercept = slopeOut = EROut = np.nan
    else:
        poIntercept = x_y_intercept(cdo, po, confinementHeight, widthT/2, True, centerPointHeight)
        slopeOut    = confinement_slope(poIntercept, po, cdo, centerPointHeight, confinementHeight, widthT, False)
        EROut       = confinement_ratio(poIntercept, widthT)
    if (isinstance(pi, int)):
        piIntercept = slopeInn = ERInn = np.nan
    else:
        piIntercept = x_y_intercept(cdi, pi, confinementHeight, widthT/2, True, centerPointHeight)
        slopeInn    = confinement_slope(piIntercept, pi, cdi, centerPointHeight, confinementHeight, widthT, False)
        ERInn       = confinement_ratio(piIntercept, widthT)
        

    return poIntercept, piIntercept, centerPointHeight, confinementHeight, slopeOut, slopeInn, EROut, ERInn
