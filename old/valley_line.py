#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 09:58:07 2024

@author: niekcollotdescury
"""

from line_functions import multiline_check, remove_man_add_section, split_ring
from smoothing import SG_smoothing
from support import create_dir


import shapely
from shapely.geometry import LineString,Polygon, Point, MultiPoint
import centerline_width
import geopandas as gpd

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def find_offset(line, width, smoothing):
    # print('find_offset')
    
    buffer_init = int((width / 2))
    
    
    intersect_line_L = multiline_check(line.offset_curve(1  * buffer_init)) # L
    intersect_line_R = multiline_check(line.offset_curve(-1 * buffer_init)) # R
    L_length = intersect_line_L.length / line.length
    R_length = intersect_line_R.length / line.length
    poly = Polygon([*list(intersect_line_L.coords), 
                                     *list(intersect_line_R.coords)[::-1]])


    buffer = buffer_init * 8
    # buffer = buffer_init * 10
    
    check = 0 
    check1 = 0
    check2 = 0

    error = False
    
    bufferList = []
    c = 0

    # plt.figure()
    # print(intersect_line_L.intersects(intersect_line_R))
    # print(intersect_line_R)
    # plt.plot(*line.xy, label = 'CL')
    # plt.plot(*poly.exterior.xy, label = 'intersect poly')
    # plt.plot(*intersect_line_L.xy, label = 'intersect L', linestyle = '--')
    # plt.plot(*intersect_line_R.xy, label = 'intersect R', linestyle = '--')
    # plt.legend()
    # plt.show()
    # l = line.segmentize(1)
    # print('start while loop', width, smoothing, buffer, len(line.coords), len(l.coords))
    if (poly.intersects(line)) & (intersect_line_L.is_empty == False) & (intersect_line_R.is_empty == False):
        while check < 2:
            
            buffer_change = 0.90
            bufferList.append(buffer)
    
            # print(f'Buffer: {buffer}, BI {buffer_init}, W {width}, SM {smoothing}, OFS {offset_side}')
            
            buffered_line_L = multiline_check(line.offset_curve(1  * buffer)).segmentize(1)
            buffered_line_R = multiline_check(line.offset_curve(-1 * buffer)).segmentize(1)
            # print(len(buffered_line.coords))
            # plt.figure()
            # plt.plot(*line.xy, label = 'CL')
            # plt.plot(*poly.exterior.xy, label = 'intersect poly')
            if (len(buffered_line_L.coords) > 2) & (len(buffered_line_R.coords) > 2):
                
                # tm = datetime.now()
                offset_L = SG_smoothing(buffered_line_L, width*smoothing)
                offset_R = SG_smoothing(buffered_line_R, width*smoothing)
                # print('SG smoothing: ',datetime.now() - tm)
                # print(f'offset while loop2 {c}_{buffer}')

                # plt.plot(*offset_L.xy, label = f'L {buffer}_counter = {c}', color = 'b')
                # plt.plot(*buffered_line_L.xy, label = f'BLL', color = 'c')
                # plt.plot(*offset_R.xy, label = f'R {buffer}_counter = {c}', color = 'b')
                # plt.plot(*buffered_line_R.xy, label = f'BLR', color = 'c')
    
                
                # plt.plot(*poly.exterior.xy, label = 'BL_poly', color = 'coral')
                # plt.plot(*intersect_line_L.xy, label = f'BL_l', color = 'firebrick')
                # plt.plot(*intersect_line_R.xy, label = f'BL_r', color = 'r')
                
                # plt.plot(*line.xy, label = f'CL', color = 'yellow')
                # print(f'offset while loop3 {c}')
                
                if (poly.intersects(offset_L)) | (poly.intersects(offset_R)):
                    # print('intersect == True')
                    buffer_change = 1.1
                    check1 = 1
                elif (poly.intersects(offset_L) == False) & (poly.intersects(offset_R) == False) & (check1 == 1):
                    # print('intersect == False & check1 == True')
                    check2 = 1
    
                # print(f'offset while loop5 {c}')
    
            check = check1 + check2
            # print('check', check)
        
            if (buffer < buffer_init) & (check1 != 1):
                print('Error: in find offset function: buffer to small without solution')
                error = True
                check = 2
            
            buffer *= buffer_change
            
            if c > 200:
                print('Error: in find offset function: while loop to long')
                error = True
                check = 2
            c +=1
            # plt.title()
            # plt.legend()
            # plt.show()
    
    else:
        error = True
        print('river boundary not correctly created')
    # print(c)
    if error == False:
        buffer = bufferList[c-1]
        offset_L = SG_smoothing(multiline_check(line.offset_curve(1  * buffer)), width*smoothing)
        offset_R = SG_smoothing(multiline_check(line.offset_curve(-1 * buffer)), width*smoothing)
        if offset_L.intersects(offset_R):
            p = shapely.ops.nearest_points(offset_L,offset_R)[0]
            dists = [p.distance(offset_L.interpolate(0, normalized = True)), p.distance(offset_L.interpolate(1, normalized = True))]
            if dists[0] > dists[1]:
                pos = 1
            else:
                pos = 0
            
            offset_L = remove_man_add_section(offset_L, p, offset_L.interpolate(pos, normalized = True))
            offset_R = remove_man_add_section(offset_R, p, offset_R.interpolate(pos, normalized = True))
    else:
        buffer = np.nan
        offset_L = LineString()
        offset_R = LineString()
    
    return offset_L, offset_R, buffer, L_length, R_length


def centerline_from_edges(L, R, line, buffer, width, smoothing):
    error = False

    cent = []
    RLs = []
   
    # plt.figure()
    # plt.title('centerline_from_edges')
    # plt.plot(*L.xy)
    # plt.plot(*R.xy)
    # plt.plot(*line.xy)

    if L.length > R.length:
        p_line = L
        d_line = R
        
    else:
        p_line = R
        d_line = L

    #ADD first point to line
    first_pd = d_line.interpolate(0, normalized = True)
    first_pp = p_line.interpolate(0, normalized = True)        
    RL = LineString([first_pd,first_pp])
    RLs.append(RL)
    cent.append(RL.interpolate(RL.length / 2))

    segments = np.linspace(start = 0,stop = 1, num = 150)
    for i in range(len(segments)):

        
        Lp = p_line.interpolate(segments[i], normalized = True)
        Rp = shapely.ops.nearest_points(Lp, d_line)[1]
        RL = LineString([Lp, Rp])
        # print(len(RL.segmentize(1).coords))
        if len(RL.segmentize(1).coords) > 3:
            RL_short = LineString(RL.segmentize(1).coords[1:-1])    
            
            if (RL_short.intersects(p_line) == False):
                cent.append(RL.interpolate(RL.length / 2))
                RLs.append(RL)


    # ADD final point
    final_pd = d_line.interpolate(1, normalized = True)
    final_pp = p_line.interpolate(1, normalized = True)
    RL = LineString([final_pd,final_pp])
    RLs.append(RL)
    cent.append(RL.interpolate(RL.length / 2))

    
    centLine = LineString(cent)
    centLine = second_centerline(centLine, L, R, line)
    # centLine = second_centerline(centLine, L, R, line)
    if (centLine.intersects(L)) | (centLine.intersects(R)):
        error == True
    
    if error == True:
        centLine = LineString()
        RLs = []
        L = LineString()
        R = LineString()
        buffer = np.nan

    return centLine, RLs,L,R, buffer


# #### second centerline 2
def second_centerline(line, off_L, off_R, river_line):
    # tm = datetime.now()
    num_segments = 150
    segments = np.linspace(0,1,num_segments)

    cents = []
    for i in range(num_segments):
        pc = line.interpolate(segments[i], normalized = True)
        pl = shapely.ops.nearest_points(pc, off_L)[1]
        pr = shapely.ops.nearest_points(pc, off_R)[1]
        dist_cl = pc.distance(pl)
        dist_cr = pc.distance(pr)

        if dist_cl > dist_cr:
            cl_line = LineString([pc, pl])
            
            new_p = cl_line.interpolate((dist_cl - dist_cr) / 2)
        else:
            cr_line = LineString([pc, pr])
            
            new_p = cr_line.interpolate((dist_cr - dist_cl) / 2)
            
        cents.append(new_p)
    new_line = LineString(cents)
    # print('Second centerline moving:', datetime.now() - tm)
    # plt.figure()
    # plt.plot(*line.xy)
    # plt.plot(*off_L.xy)
    # plt.plot(*off_R.xy)
    # plt.plot(*river_line.xy)
    # plt.plot(*new_line.xy)
    # tm = datetime.now()
    new_line = SG_smoothing(new_line, 200, True)
    # print('Second centerline smoothing:', datetime.now() - tm)
    # plt.plot(*new_line.xy)
    # plt.show()
    return new_line



# #### Call all functions (River_dir) 

# one time
def river_dir(df,df_node, smoothing):
    def no_returns():
        return (shapely.geometry.LineString(), [], shapely.geometry.LineString(), shapely.geometry.LineString(), np.nan)
        
    width = int(df_node.max_width.mean())
    
    # tm = datetime.now()
    offset_L, offset_R,buffer, L_length, R_length, = find_offset(df.geometry, width,smoothing)
    # print('Find offset time:', datetime.now() - tm)
    # print(buffer)
    # offset_R,buffer_R = find_offset(df.geometry, width,smoothing, -1, buffer_L)

    # plt.plot(*offset_L.xy, label = 'L', color = 'r')
    # plt.plot(*offset_R.xy, label = 'R', color = 'r')
    # plt.plot(*df.geometry.xy, label = 'center', linestyle = '--', color = 'b')
    # plt.legend()
    # plt.show()
    
    # tm = datetime.now()
    if np.isnan(buffer) == False:
        dir_centerline, RLs,offset_L,offset_R, buffer = centerline_from_edges(offset_L, offset_R, df.geometry,buffer, width, smoothing)
    else:
        dir_centerline ,RLs, offset_L, offset_R, buffer = no_returns()
    # print('Centerline time:', datetime.now() - tm)
    return (dir_centerline, offset_L, offset_R, RLs, buffer, L_length, R_length)



# Inflection point method ######################################################
def centerline_width_apply(L, R, id,projection, projection_ee, riverLine, directory):
    gdf_L = gpd.GeoDataFrame({'geometry': [L]}, crs = projection)
    gdf_L = gdf_L.to_crs(projection_ee).iloc[0].geometry

    gdf_R = gpd.GeoDataFrame({'geometry': [R]}, crs = projection)
    gdf_R = gdf_R.to_crs(projection_ee).iloc[0].geometry

    create_dir(directory + 'centerline_temp_csv')
    
    # temp_df = pd.DataFrame({'llat':llat,'llon':llon,'rlat':rlat,'rlon':rlon})
    temp_df = pd.DataFrame()
    
    # Note this column has 4 rows of values:
    slicing = [0,415]
    pdl = pd.DataFrame({'llat': gdf_L.xy[1][slicing[0]:slicing[1]], 'llon': gdf_L.xy[0][slicing[0]:slicing[1]]})
    pdr = pd.DataFrame({'rlat': gdf_R.xy[1][slicing[0]:slicing[1]], 'rlon': gdf_R.xy[0][slicing[0]:slicing[1]]})

    pdl = pd.DataFrame({'llat': gdf_L.xy[1], 'llon': gdf_L.xy[0]})
    pdr = pd.DataFrame({'rlat': gdf_R.xy[1], 'rlon': gdf_R.xy[0]})
    temp_df = pd.concat([temp_df, pdl], axis=1) 
    temp_df = pd.concat([temp_df, pdr], axis=1) 
    
    temp_df.to_csv(directory + f'centerline_temp_csv/{id}_temp.csv', sep=',', index = False)

    river_object = centerline_width.riverCenterline(csv_data=directory + f'centerline_temp_csv/{id}_temp.csv', 
                                                    ellipsoid = 'WGS84', interpolate_data=True)

    
    centerline = LineString(river_object.centerlineSmoothed)
    
    df_pee = gpd.GeoDataFrame({'geometry': [centerline]}, crs = projection_ee)
    centerline = df_pee.to_crs(projection).iloc[0].geometry
    return centerline

    

def inflection_valley(row, apex_points, smoothing, id,projection, projection_ee, directory, plot = False):
    line = row.geometry
    width = row.node_mwm
    L = line.buffer(width / 2)

    apex_points.insert(0, Point(line.coords[0]))
    apex_points.append(Point(line.coords[-1]))
    
    hull = shapely.concave_hull(MultiPoint(apex_points))
    
    
    endLineSeparator = shapely.GeometryCollection([apex_points[0].buffer(1), apex_points[-1].buffer(1)])
    
    hullSplit = split_ring(hull.exterior, endLineSeparator)
    
    if len(hullSplit.geoms) > 2:
        print('check number of splits and take biggest two')
    else:
        side1 = hullSplit.geoms[0]
        side2 = hullSplit.geoms[1]


    if Point(side1.coords[0]).distance(Point(side2.coords[0])) > 1e3:
        side2 = shapely.reverse(side2)
        offset_side = -1
    else: offset_side = 1
    
    smoothing_window = smoothing*width               # define smoothing window based on mean node max width
    length_ratio     = line.length / row.min_len     # Length ratio reach length and minimum required length
    if length_ratio < 1.0: # If length ratio is below 1 the smoothing window is multiplied by the ratio
        smoothing_window *= length_ratio

    s1Inter = SG_smoothing(side1, smoothing_window).intersects(L)
    s2Inter = SG_smoothing(side2, smoothing_window).intersects(L)

    
    buffer_width = width / 2
    
    while (s1Inter) | (s2Inter):

        s1 = side1.offset_curve(buffer_width              , join_style = 'round')
        s2 = side2.offset_curve(buffer_width * offset_side, join_style = 'round')
      
        
        
        buffer_width *= 1.1

        
        side1Smooth = SG_smoothing(s1, smoothing_window)
        side2Smooth = SG_smoothing(s2, smoothing_window)
        
        s1Inter = side1Smooth.intersects(L)
        s2Inter = side2Smooth.intersects(L)

    # dir_centerline, RLs,offset_L,offset_R, buffer = centerline_from_edges(side1Smooth, side2Smooth, line ,buffer_width, width, smoothing)
    dir_centerline = centerline_width_apply(side1Smooth, side2Smooth, id,projection, projection_ee, line , directory)

    sinuosity = line.length / dir_centerline.length
    if plot == True:
        plt.figure()
        plt.plot(*L.exterior.xy, 'g')
        # plt.plot(*hull.exterior.xy)
        plt.plot(*line.xy, 'orange')
        plt.plot(*side1.xy, 'yellow')
        plt.plot(*side2.xy, 'yellow')
        plt.plot(*dir_centerline.xy, 'b')
    
        plt.plot(*side1Smooth.xy, 'r')
        plt.plot(*side2Smooth.xy, 'r')
        sinus_title = f'Sinuosity: {np.round(sinuosity, 3)}'
        lineInfo_title = f'Length: {int(line.length)}, Width: {int(width)}'
        plt.title(f'{id}\n{lineInfo_title}\n{sinus_title}')
        plt.axis('equal')
        plt.axis('off')
        plt.show()

    return dir_centerline
# id = 31101000031

# D = df[df.reach_id == id]
# R = D.iloc[0]

# sin,reg, new_df, infs, apex, apex_points = inflection_points(id, D, False, 20,10 , True)
# inflection_valley(R.geometry,apex_points, R.node_mwm, 50, R.reach_id, True)

