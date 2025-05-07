#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 09:58:06 2024

@author: niekcollotdescury
"""
# import packages
import shapely
import numpy as np
import math as m
import geopandas as gpd

# import partial packages
from shapely.geometry import LineString, Point, Polygon, MultiLineString, MultiPolygon

def multiline_check(line):
    if line.geom_type == 'MultiLineString':
        line = shapely.ops.linemerge(line)
    
    if line.geom_type == 'MultiLineString':
        l = np.zeros(len(line.geoms))
        for i in range(len(line.geoms)):
            l[i] = line.geoms[i].length
        line = line.geoms[l.argmax()]
    
    return line

def remove_man_add_section(line, node, endP):
        
    P = line.interpolate(line.project(node))

    splits = shapely.ops.split(line, P.buffer(1))

    dist = np.zeros(len(splits.geoms))
    lens = np.zeros(len(splits.geoms))
    for s in range(len(splits.geoms)):
        dist[s] = endP.distance(splits.geoms[s])
        lens[s] = splits.geoms[s].length
        

    
    distSort = dist.argsort()
    if lens.argmin() == dist.argmax():
        distIndex = distSort[-2]
    else:
        distIndex = distSort[-1]
    
    newLine = splits.geoms[distIndex]

    return newLine

def getExtrapoledLine(p1,p2, w, wm, single = False, factor_distance = 'Factor/'):
    """Creates a line extrapoled in p1->p2 direction\n
    input: 
    - p1: starting Point, tuple or list of coordinate pair.
    - p2: end point, tuple or list of coordinate pair.
    - w:  ratio or extension of linestring for distance from p1
    - wm: 2 or 1. Only used if factor is "Factor"
        - 2: linestring is increased with factor w or distance.
        - 1: linestring is increased with factor w/2.
    - single: extend line in direction of p2 or both directions from p1.
    - factor_distance: extend line based on factor of w/wm or use w as distance (Factor or Distance)\n
    return:\n
        extended linestring
    """
    if factor_distance == 'Factor':
        EXTRAPOL_RATIO = w*(wm/2)  # old method
    elif factor_distance == 'Distance':
        dist = Point(p1).distance(Point(p2))
        EXTRAPOL_RATIO = w / dist
    else:
        print('ERROR give correct value for factor_distance')
    
    if single == True:
        a = p1
        b = (p1[0]+EXTRAPOL_RATIO*(p2[0]-p1[0]), p1[1]+EXTRAPOL_RATIO*(p2[1]-p1[1]) )
        # b = (p2[0]+EXTRAPOL_RATIO*(p1[0]-p2[0]), p2[1]+EXTRAPOL_RATIO*(p1[1]-p2[1]) )
    else:
        a = (p1[0]+EXTRAPOL_RATIO*(p2[0]-p1[0]), p1[1]+EXTRAPOL_RATIO*(p2[1]-p1[1]) )
        b = (p2[0]+EXTRAPOL_RATIO*(p1[0]-p2[0]), p2[1]+EXTRAPOL_RATIO*(p1[1]-p2[1]) )
    return LineString([a,b])

def azimuth_coords(x1,y1,x2,y2):
    '''azimuth between 2 shapely points (interval 0 - 360)'''
    angle = np.arctan2(x2 - x1, y2 - y1)
    return np.degrees(angle) if angle >= 0 else np.degrees(angle) + 360

def azimuth(point1, point2):
    '''azimuth (orientation) between 2 shapely points (interval 0 - 360)\n
    input:
    - point1/point2: Shapely Point\n
    output: orientation in degrees
    - 
    '''
    angle = np.arctan2(point2.x - point1.x, point2.y - point1.y)
    return np.degrees(angle) if angle >= 0 else np.degrees(angle) + 360

def split_ring(ring, split):
    """Split a linear ring geometry, returns a [Multi]LineString

    See my PostGIS function on scigen named ST_SplitRing
    """
    valid_types = ('MultiLineString', 'LineString', 'GeometryCollection')
    if not hasattr(ring, 'geom_type'):
        raise ValueError('expected ring as a geometry')
    elif not hasattr(split, 'geom_type'):
        raise ValueError('expected split as a geometry')
    if ring.geom_type == 'LinearRing':
        ring = LineString(ring)
    if ring.geom_type != 'LineString':
        raise ValueError(
            'ring is not a LinearRing or LineString, found '
            + str(ring.geom_type))
    elif not ring.is_closed:
        raise ValueError('ring is not closed')
    elif split.is_empty:
        return ring
    elif not split.intersects(ring):
        # split does not intersect ring
        return ring
    if split.geom_type == 'LinearRing':
        split = LineString(split)
    if split.geom_type not in valid_types:
        raise ValueError(
            'split is not a LineString-like or GeometryCollection geometry, '
            'found ' + str(split.geom_type))

    intersections = ring.intersection(split)
    if intersections.is_empty:
        # no intersections, returning same ring
        return ring
    elif intersections.geom_type == 'Point':
        # Simple case, where there is only one line intersecting the ring
        result = Polygon(ring).difference(split).exterior
        # If it is a coordinate of the ring, then the ring needs to be rotated
        coords = result.coords[:-1]
        found_i = 0
        for i, c in enumerate(coords):
            if Point(c).almost_equals(intersections):
                found_i = i
                break
        if found_i > 0:
            result = Polygon(coords[i:] + coords[:i]).exterior
        if result.interpolate(0).distance(intersections) > 0:
            raise Exception(
                'result start point %s to intersection %s is %s' %
                (result.interpolate(0), intersections,
                 result.distance(intersections)))
        elif result.geom_type != 'LinearRing':
            raise Exception(
                'result is not a LinearRing, found ' + result.geom_type)
        elif not result.is_closed:
            raise Exception('result is not closed')
        return LineString(result)

    difference = ring.difference(split)
    if difference.geom_type != 'MultiLineString':
        raise ValueError(
            'expected MultiLineString difference, found '
            + difference.geom_type)

    start_point = ring.interpolate(0)
    if start_point.distance(intersections) == 0:
        # special case: start point is the same as an intersection
        return difference

    # Otherwise the line where the close meets needs to be fused
    fuse = []
    parts = list(difference.geoms)
    for ipart, part in enumerate(parts):
        if part.intersects(start_point):
            fuse.append(ipart)
    if len(fuse) != 2:
        raise ValueError('expected 2 geometries, found ' + str(len(fuse)))
    # glue the last to the first
    popped_part = parts.pop(fuse[1])
    parts[fuse[0]] = shapely.ops.linemerge([parts[fuse[0]], popped_part])
    return MultiLineString(parts)

def remove_self_intersect(Line):
    if shapely.is_simple(Line) == False:
        mls = shapely.ops.unary_union(Line)

        coords_list = []
        for l in mls.geoms:
            if len(shapely.ops.polygonize(l)) ==0:
                coords_list.extend(l.coords)
        Line = shapely.LineString(coords_list)
    return Line

def create_angled_lines(line, dist = 50,max_angle = 10, return_angle = False):
    """
    Create line that deviates with set distance from main line.\n
    input:
    - line: original line
    - dist: distance of the opposite side
    - return_angle: boolean to return angle\n
    return:
    - line for positive angle and negative angle
    - angle: optional return of the angle between original line and created lines
    """
    lineLen = line.length


    alpha = np.rad2deg(np.arcsin((dist)/(lineLen)))
    while alpha > max_angle:
        dist *= 0.9
        alpha = np.rad2deg(np.arcsin((dist)/(lineLen)))


    # get first and last coords of original line
    p1 = line.coords[0]
    p2 = line.coords[-1]

    # get the angle of the original line to a zero degree line
    lineAngle = m.atan2(p2[1] - p1[1], p2[0]  - p1[0])
    angle     = m.radians(alpha)
    

    # calculate position of end points and create line. Repeat for negative angle
    endxp = p1[0] + lineLen * m.cos(angle + lineAngle)
    endyp = p1[1] + lineLen * m.sin(angle + lineAngle)
    angled_linep = LineString([p1, (endxp, endyp)])
    

    endxm = p1[0] + lineLen * m.cos((angle *-1) + lineAngle)
    endym = p1[1] + lineLen * m.sin((angle *-1) + lineAngle)
    angled_linem = LineString([p1, (endxm, endym)])
    if return_angle == True:
        return (angled_linep, angled_linem, alpha)
    else:
        return (angled_linep, angled_linem)

def get_points_along_linestring(line, spacing = 20, separate_xy = False):
    num_points = int(line.length / spacing) + 1
    # Compute the total length of the LineString
    total_length = line.length
    
    # Generate distances at which points will be located
    distances = np.linspace(0, total_length, num_points)
    
    # Interpolate points at the specified distances
    points = [list(line.interpolate(distance).coords[0]) for distance in distances]
    return points

def angle_diff(L : "LineString", stepSize :'int' = 50, plot : 'bool' = False):
    '''
    Function that calculates the difference in angle between concecutive points. 
    Number of points determined by the step input value. 
    input:
        L: River Centerline 
        total_length: total length of river segment
        step: Percent stepsize of between concecutive points 
    '''
    points = get_points_along_linestring(L,50)

    azs = []

    for i in range(len(points) -2):
        P_m = points[i]
        P_t = points[i+1]
        P_p = points[i+2]
        # Calculate direction between 
        fwd_azimuth_minus = azimuth_coords(P_m[0], P_m[1], P_t[0], P_t[1])
        fwd_azimuth_plus  = azimuth_coords(P_t[0], P_t[1], P_p[0], P_p[1])
        az_diff = abs(fwd_azimuth_minus - fwd_azimuth_plus)
        if az_diff > 180:
            az_diff = abs(az_diff - 360)
            
        
        azs.append(az_diff)
        
        
    # remove nan values
    return np.nanmean(azs)

def curvature(L,coord_selection = True, average = True):
    """
    Calculate Curvature of LineString based on first and second derivative
    input:
        L: LineString
        plot: Boolean for plotting (default = False)
    output:
        averaged value for the curvature 
    """
    # Select coordinates for linestring. Either use exisitng coordinates or select coordinates
    # based on spacing
    if coord_selection == True:
        points = get_points_along_linestring(L,10)
        if len(points) == 1:
            points = get_points_along_linestring(L,1)
        if len(points) == 1:
            return np.nan
    else:
        points = np.array(L.coords)
        if len(points) <= 1:
            L = L.segmentize(1)
        if len(points) <= 1:
            return np.nan
        points = np.array(L.coords)

    ###########################
    # Calculate Curvature
    x = [p[0] for p in points]
    y = [p[1] for p in points]

    # Calculate first derivatives using finite differences
    dx = np.gradient(x)
    dy = np.gradient(y)

    # Calculate second derivatives using finite differences
    ddx = np.gradient(dx)
    ddy = np.gradient(dy)

    # Calculate curvature using the formula
    curvature = (dx * ddy - dy * ddx) / (dx**2 + dy**2)**(3/2)

    # Average curvature for Line segment
    if average == False:
        curve = curvature
    else:
        curve = np.abs(np.nanmean(curvature))
    return curve

# def get_bend_width(line: 'LineString', DFN: 'gpd.GeoDataFrame'):
    """
    Get the width of individual bends determined by the proximity of nodes to the bend
    input: 
        line: linestring 
        reaches: reach_ids associated with the linestring
        DFN: node dataframe
        crs: local crs
        plot: boolean for plotting (default = True)
    return:
        bendwidth in meters
    """
    reach_nodes = DFN.copy()
    dfBend = gpd.GeoDataFrame({'ID':[1], 'geometry':line.buffer(200)}, crs = reach_nodes.crs)

    lineNodeJoin = gpd.sjoin_nearest(dfBend, reach_nodes)
    bendNodes    = reach_nodes[reach_nodes.node_id.isin(lineNodeJoin.node_id.values)]
    bendWidthMax = bendNodes.max_width.mean()
    bendWidth    = bendNodes.width.mean()
    
    return bendWidth, bendWidthMax

def get_bend_width(line, bendLine, dfN, dfR):
    
    # region Fold Start
    start = line.project(Point(bendLine.coords[0]))
    end   = line.project(Point(bendLine.coords[-1]))
    
    nodes = dfN[(dfN['linePos'] > start) & (dfN['linePos'] < end)]
    if nodes.shape[0] == 0:
        nodeDistance = dfN['linePos'] - ((end-start) / 2)
        nodes        = dfN.iloc[nodeDistance.argmin()]

    # endregion Fold End
    meanWidth    = nodes['width'].mean()
    meanMaxWidth = nodes['max_width'].mean()
    
    if (np.isnan(meanWidth)) | (np.isnan(meanMaxWidth)):
        if isinstance(nodes['reach_id'], np.int64):
            nodeReaches = [nodes['reach_id']]
        else:
            nodeReaches = nodes['reach_id']
        reaches = dfR[dfR['reach_id'].isin(nodeReaches)]
        if (np.isnan(meanWidth)):
            meanWidth = reaches['width'].mean()
        if np.isnan(meanMaxWidth):
            meanMaxWidth = reaches['max_width'].mean()

    return meanWidth, meanMaxWidth

def get_bend_dist_out(line, bendLine, dfN):
    # region Fold Start
    start = line.project(Point(bendLine.coords[0]))
    end   = line.project(Point(bendLine.coords[-1]))
    
    nodes = dfN[(dfN['linePos'] > start) & (dfN['linePos'] < end)]
    if nodes.shape[0] == 0:
        nodeDistance = dfN['linePos'] - ((end-start) / 2)
        nodes        = dfN.iloc[nodeDistance.argmin()]

    # endregion Fold End
    bendDistOut    = nodes['dist_out'].max()

    return bendDistOut


# def coord_distances(line):
    coords = np.array(line.coords)

    # Compute segment lengths
    diffs              = np.diff(coords, axis=0)  # Vector differences between points
    segment_lengths    = np.linalg.norm(diffs, axis=1)  # Euclidean distances
    cumulative_lengths = np.insert(np.cumsum(segment_lengths), 0, 0)  # Include start at 0

    coordDist = np.zeros(len(coords))
    for i,coord in enumerate(coords):
    # Find index of target point in coordinates
        idx = np.where((coords == (coord[0], coord[1])).all(axis=1))[0]

        if idx.size > 0:
            coordDist[i] = cumulative_lengths[idx[0]]  # Distance along the LineString

        else:
            print("Point not found in LineString.")
    assert len(coords) == len(coordDist), 'positions of coords not equal in length to coords'
    return coordDist

# def create_bend_line(P1, P2, coordDist, line):

    P1Ind = int(abs(coordDist - line.project(P1)).argmin())
    P2Ind = int(abs(coordDist - line.project(P2)).argmin())

    lineCoords = line.coords

    newLine = LineString(lineCoords[P1Ind:P2Ind+1])

    return newLine

def create_bend_line(inf_line, riv_line):
    """
    create lines corresponding to inflection sections
    input:
        inf_line: inflection line sections
        riv_line: river centerline
    output:
        centerline associated with inflection line sections
    """

    # cut the river line in sections corresponding to the inflection section
    p1 = inf_line.interpolate(inf_line.project(Point(inf_line.coords[0 ])))
    p2 = inf_line.interpolate(inf_line.project(Point(inf_line.coords[-1])))
    p1B = p1.buffer(10)
    p2B = p2.buffer(10)

    p = MultiPolygon([p1B,p2B])
    vectorSplit = shapely.ops.split(riv_line, p)

    #add assert for number of splits??????????

    # loop over the split sections and find the sections with the minimum distance to the start end point of the sections
    splitDist = np.zeros(len(vectorSplit.geoms))
    for j in range(len(vectorSplit.geoms)):
        splitLine = vectorSplit.geoms[j]
        splitDist[j] = splitLine.distance(p1) + splitLine.distance(p2)
    alongSectionId = splitDist.argmin()

    # extract correct split section
    alongSection = vectorSplit.geoms[alongSectionId]
    L = LineString(list([p1]) + list(alongSection.coords) + list(p2.coords))

    return L 


def checkIntersection(line1, line2):
    line1 = LineString([line1.interpolate(0.1,normalized = True),
                        line1.interpolate(1  ,normalized = True)])

    if line1.intersects(line2):
        p1 = line1.intersection(line2) 

        if (p1.geom_type == 'MultiPoint'):
            p1Sorted = sorted(p1.geoms, key=lambda seg: line1.project(Point(seg.coords[0])))
            p1 = p1Sorted[0]
        elif (p1.geom_type == 'LineString'):
            p1 = Point(p1.coords[0])
        elif (p1.geom_type == 'MultiLineString'):
            p1 = Point(p1.geoms[0].coords[0])
    else:
        p1 = np.nan

    return p1

def adjust_confinement_line(line, centerline):
    interPoint = checkIntersection(line, centerline)

    if isinstance(interPoint, Point):
        interPoint  = line.interpolate(line.project(interPoint))
        vectorSplit = shapely.ops.split(line, interPoint.buffer(10))
        vectorSplit = sorted(vectorSplit.geoms, key=lambda seg: line.project(Point(seg.coords[0])))
        line = vectorSplit[0]
    return line

def localCRSReachLength(DF, col, projection):
    df = DF.copy()
    for CRS in df.localCRS.unique():
        rows = df[df.localCRS == CRS]
        rows = rows.to_crs(CRS)
        df.loc[rows.index, col] = rows.geometry.length
    df = df.to_crs(projection)
    return df