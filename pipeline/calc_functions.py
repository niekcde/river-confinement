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

def confinement_slope(intercept, profile, distance,centerPointHeight, confHeight, width, normalize = True):
    heightDiff = confHeight - centerPointHeight
    maxSlope = heightDiff / (width / 2)

    if np.isnan(intercept):
        if np.max(distance) < (width / 2):
            return np.nan
        else:
            slope = slope_with_intercept(distance, profile, centerPointHeight)
    else:
        
        
        if intercept < (width / 2):
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
        if intercept < (width/2):
            ratio = 0.5
    return ratio


def _clean_profile_pair(profile, distance):
    if isinstance(profile, int):
        return np.array([], dtype=float), np.array([], dtype=float)
    if isinstance(distance, int):
        return np.array([], dtype=float), np.array([], dtype=float)

    try:
        profile_arr = np.asarray(profile, dtype=float)
        distance_arr = np.asarray(distance, dtype=float)
    except Exception:
        return np.array([], dtype=float), np.array([], dtype=float)

    if profile_arr.size == 0 or distance_arr.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float)

    valid = (
        np.isfinite(profile_arr)
        & np.isfinite(distance_arr)
        & (profile_arr != 99999)
        & (distance_arr != 99999)
    )
    return profile_arr[valid], distance_arr[valid]

def get_low_center_point(po, pi, cdo, cdi, width):
    po, cdo = _clean_profile_pair(po, cdo)
    pi, cdi = _clean_profile_pair(pi, cdi)

    poRiv = po[cdo <= (width / 2)] if po.size else np.array([], dtype=float)
    piRiv = pi[cdi <= (width / 2)] if pi.size else np.array([], dtype=float)

    mins = []
    if poRiv.size:
        mins.append(np.min(poRiv))
    if piRiv.size:
        mins.append(np.min(piRiv))

    if len(mins) == 0:
        return np.nan
    return np.min(mins)

# Function to convert mm to inches
def mm_to_inch(mm):
    return mm / 25.4

def log_var(df, cols):
    """ Apply log10 to columns, with an error adjustment. 
    Error adjustment calculated as min values divided by 10.\n
    input:
    - df: input dataset
    - cols: columns to be logged\n
    output:
    - dataset
    """
    for c in cols:
        eps = df[df[c] >0][c].min()/10
        df[f'{c}Log'] = np.log10(df[c] + eps)
    return df
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

def scale(df, scaleCols):
    """
    Apply standard scaler to columns in scaleCols list.
    All NA values dropped from selected columns
    return scaled columns in array format.
    """
    df = df[scaleCols].copy()
    df = df.dropna()
    X  = df[scaleCols]
    # === Standardize features ===
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    return X_scaled

def stratified_sample(X, labels, min_labels, n, rs = 42):
    """
    Take stratified sample from input dataset\n
    Input:
    - X: input data to be stratified (arrays or pandas data series)
    - labels: labels corresponding to data
    - min_labels: minimal number of unique labels
    - n: number of samples to taken from dataset
    - rs: random_state assigned for sampling (default 42)
    """
    # test if desired number of labels is present and not all labels are missing
    if len(np.unique(labels)) < min_labels or all(labels == -1):
        print('sampling not possible due to mismatch unique samples\
            and minimal number of samples')
        return None, None
    
    X_sample, _, y_sample, _ = train_test_split(X, labels, train_size=n, 
                                                stratify=labels, random_state=rs)
    return X_sample, y_sample

def confinement_values(po, pi, cdo, cdi,widthW, widthT, factor):
    """
    Calculate confinement slope and ratio of intercept with river width.\n
    input:
    - po/pi: elevation profile for outer and inner bend
    - cdo/cdi: distance from river centerline for outer and inner bend
    - widthW: wetted river width
    - widthM: Total river width
    - factor: factor that multiplied with the river width determines the confinement height\n
    output:
    - po/pi Intercept: intercept distance for outer and inner bend
    - centerPointHeight: height of the centerpoint of the river determmined by the lowest point within river bounds
    - confinementHeight: valley edge height based on factor multiplied with river width above river centerpoint
    - slope Out/Inn: normalized outer and inner slope
    - ER Out/Inn: Ratio of valley edge intercept and river width for outer and inner bend
    """

    po_clean, cdo_clean = _clean_profile_pair(po, cdo)
    pi_clean, cdi_clean = _clean_profile_pair(pi, cdi)

    centerPointHeight = get_low_center_point(po_clean, pi_clean, cdo_clean, cdi_clean, widthW)
    if np.isnan(centerPointHeight):
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    confinementHeight = (widthW * factor) + centerPointHeight

    if po_clean.size == 0:
        poIntercept = slopeOut = EROut = np.nan
    else:
        poIntercept = x_y_intercept(cdo_clean, po_clean, confinementHeight, widthW/2, True, centerPointHeight)
        slopeOut    = confinement_slope(poIntercept, po_clean, cdo_clean, centerPointHeight, confinementHeight, widthW, False)
        EROut       = confinement_ratio(poIntercept, widthW)
    if pi_clean.size == 0:
        piIntercept = slopeInn = ERInn = np.nan
    else:
        piIntercept = x_y_intercept(cdi_clean, pi_clean, confinementHeight, widthW/2, True, centerPointHeight)
        slopeInn    = confinement_slope(piIntercept, pi_clean, cdi_clean, centerPointHeight, confinementHeight, widthW, False)
        ERInn       = confinement_ratio(piIntercept, widthW)
        

    return poIntercept, piIntercept, centerPointHeight, confinementHeight, slopeOut, slopeInn, EROut, ERInn
