#!/usr/bin/env python
# coding: utf-8

# # Code

# In[1]:


# Valley slope --> ?
# Slope --> CHECK
# Sinuosity -- > CHECK??
# directional value?
#Confined unconfined? --> ration between channel width and valley bottom width
# width ratio --> CHECK
# number of islands --> number of channels??


# In[2]:


# Import Package 
import numpy as np


import io
import requests
import ee
import rioxarray


import re
import itertools 

from rasterio.mask import mask

from shapely import geometry, ops
from shapely.geometry import LineString, Point, Polygon
import shapely
from shapely.ops import nearest_points
import glob

# import Read_inputs as ReadInputs
# import river_utils as ru

import scipy as sc
from scipy import misc
# from scipy.ndimage import label, binary_dilation, generate_binary_structure

import matplotlib.pyplot as plt

import geopandas as gpd
import pandas as pd

import rasterio
from rasterio import features

from pylab import *
import importlib

from osgeo import gdal

from tqdm import tqdm

import os
# from windrose import WindroseAxes
from datetime import datetime

import math as m

from shapelysmooth import taubin_smooth

import geemap
from pandas import read_csv

# from sklearn.linear_model import LinearRegression
# from sklearn.linear_model import LinearRegression
# from sklearn.preprocessing import PolynomialFeatures
# from sklearn.metrics import mean_squared_error, r2_score
from matplotlib import pyplot
import numpy


import ee

max_gdal_cache_gb=64
gdal.SetCacheMax(int(max_gdal_cache_gb * 1e9))


# conda install jupyter
# conda install shapely
# conda install geopandas


# In[3]:


directory = '/Users/niekcollotdescury/Desktop/PhD/RiverWidth/'


# # __Width method__

# In[3]:


def width_ratio(DF_vector, DF_node):
    ##### REMOVE MANUALLY ADDED VALS FOR NODE
    df_vector = DF_vector.copy()
    df_node = DF_node.copy()

    df_node = df_node[df_node.manual_add == 0]
    
    df_node['widthRatio'] = df_node.width / df_node.max_width
    df_node = df_node[df_node.reach_id.isin(df_vector.reach_id.values)]
    df = df_node.groupby('reach_id', as_index = False)['widthRatio'].mean()

    df_v = df_vector.merge(df, how = 'left', on = 'reach_id')
    
    return df_v


# #### Testing

# In[50]:


df = width_ratio(df_new_segments, df_node)
df

# pd.DataFrame.groupby?


# # __Inflection points & regularity (sinuosity)__

# ## OUD

# In[5]:


def inflection_points(ids, df_vector, plot = True, s = 10, degree = 1, end_points = False):
    #### ADJUST WIDTHHHHHHH   
    row = df_vector[df_vector.reach_id == ids].iloc[0]
    
    #smooth the line segment: factor* the segment width. (creates a line width coordinates every meter)
    vector = SG_smoothing(row.geometry, s*int(row.node_mwm))
    vec_coords       = vector.coords # coordinates in the line segment

    # empty lists to be used in loop
    crosses          = []
    inflections      = []
    inflections_true = []
    inf_loc          = []
    fwd              = []

    # Counter for the inflection points
    counter          = 0

    # loop over all the coordinates minus two. the minus two is to be able to find coordinates before and after the selected point
    for c in range(len(vec_coords) -2):
        
        # Select the target coordinate pair and select point before and after
        target_min  = vec_coords[c]
        target      = vec_coords[c+1]
        target_plus = vec_coords[c+2]
        
        # Calcilate direction between 
        fwd_azimuth_minus = azimuth_coords(target_min[0], target_min[1],target[0], target[1])
        fwd_azimuth_plus = azimuth_coords(target[0], target[1],target_plus[0], target_plus[1])
        fwd.append((fwd_azimuth_minus + fwd_azimuth_plus) / 2)
        
        v1 = np.array(target) - np.array(target_min)
        v2 = np.array(target_plus) - np.array(target)
        crosv = (v1[0] * v2[1]) - (v1[1] * v2[0])

        crosses.append(crosv)
        
        if np.sign(crosses[c]) != np.sign(crosses[c -1]):
            inflections.append(Point(target))
            
            inf_loc.append(c)
            if len(inflections) > 1:
                fwd_diff = np.abs(np.abs(fwd[inf_loc[counter]]) - np.abs(fwd[inf_loc[counter-1]]))
                # print(fwd_diff)
                if fwd_diff > degree:
                    inflections_true.append(Point(target))
            else:
                inflections_true.append(Point(target))
            counter +=1

    start = Point(vector.coords[0])
    end   = Point(vector.coords[-1])
    if end_points == True:
        inflections_true.insert(0, start)
        inflections_true.append(end)

    if len(inflections_true) > 1:
        inf_line = LineString(inflections_true)

        vector_length = vector.project(inflections_true[-1])
        inflection_length = inf_line.length
        if (vector_length == 0) | (inflection_length == 0):
            sinuosity = 0
        else:
            sinuosity = (vector_length / inflection_length) - 1
    
        sin = [] 
        for i in range(len(inflections_true) - 1):
    
            p1 = inflections_true[i]
            p2 = inflections_true[i+1]
            along_dist = vector.project(p2)-vector.project(p1)
            straight_dist = LineString([p1,p2]).length
            sin.append(along_dist / straight_dist)
        
        sin_reg = np.std(sin)
        sin_mean = np.mean(sin)
        # print(sin_mean, sin_reg)
    else:
        sinuosity = 0.00
        sin_reg   = 0.00

    # plt.show()
    if plot == True:
        plt.figure(figsize=  [6,4])
        plt.title(f'id: {ids}\nSinuosity: {np.round(sinuosity,8)}_{np.round(sin_reg,8)}\nWindow: {s}, width: {row.node_mwm}')
        plt.plot(*row.geometry.xy)
        plt.plot(*vector.xy)
        
        if len(inflections_true) > 1:
            plt.plot(*inf_line.xy)
        for p in inflections_true:
            plt.scatter(*p.xy)
        plt.axis('equal')
        plt.axis('off')
        plt.show()

    return(sinuosity, sin_reg)
    # return(sinuosity, sin_reg, inf_line)


# ## NEW

# In[5]:


def inflection_points(ids, df_vector, plot = True, s = 10, degree = 1, end_points = False):  
    row = df_vector[df_vector.reach_id == ids].iloc[0]
    
    #smooth the line segment: factor* the segment width. (creates a line width coordinates every meter)
    vector     = SG_smoothing(row.geometry, s*int(row.node_mwm), False)
    vector     = vector.simplify(0.5, preserve_topology=True)
    vec_coords = vector.coords # coordinates in the line segment

    # empty lists to be used in loop
    crosses          = []
    inflections      = []
    inflections_true = []
    inf_loc          = []
    fwd              = []

    # Counter for the inflection points
    counter          = 0

    tm = []
    t = []
    tp = []
    # loop over all the coordinates minus two. the minus two is to be able to find coordinates before and after the selected point
    for c in range(len(vec_coords) -2):
        
        # Select the target coordinate pair and select point before and after
        target_min  = vec_coords[c]
        target      = vec_coords[c+1]
        target_plus = vec_coords[c+2]
        
        # Calcilate direction between 
        fwd_azimuth_minus = azimuth_coords(target_min[0], target_min[1],target[0], target[1])
        fwd_azimuth_plus  = azimuth_coords(target[0], target[1],target_plus[0], target_plus[1])
        fwd.append((fwd_azimuth_minus + fwd_azimuth_plus) / 2)
        
        v1 = np.array(target) - np.array(target_min)
        v2 = np.array(target_plus) - np.array(target)
        crosv = (v1[0] * v2[1]) - (v1[1] * v2[0])

        crosses.append(crosv)
        
        tm.append(np.array(target_min))
        t.append(np.array(target))
        tp.append(np.array(target_plus))
        
        if np.sign(crosses[c]) != np.sign(crosses[c -1]):
            inflections.append(Point(target))
            
            inf_loc.append(c)
            if len(inflections) > 1:
                fwd_diff = np.abs(np.abs(fwd[inf_loc[counter]]) - np.abs(fwd[inf_loc[counter-1]]))
                # print(fwd_diff)
                if fwd_diff > degree:
                    inflections_true.append(Point(target))
            else:
                inflections_true.append(Point(target))
        
            counter +=1
    new_df = pd.DataFrame({'id':np.arange(0,len(vec_coords)-2), 'Cross':crosses,
                          'target_min': tm, 'target':t, 'target_plus':tp})
    new_df.loc[:,'v1'] = new_df['target'] - new_df['target_min']
    new_df.loc[:,'v2'] = new_df['target_plus'] - new_df['target']
    start = Point(vector.coords[0])
    end   = Point(vector.coords[-1])
    if end_points == True:
        inflections_true.insert(0, start)
        inflections_true.append(end)

    if len(inflections_true) > 1:
        inf_line = LineString(inflections_true)

        vector_length = vector.project(inflections_true[-1])
        inflection_length = inf_line.length
        if (vector_length == 0) | (inflection_length == 0):
            sinuosity = 0
        else:
            sinuosity = (vector_length / inflection_length) - 1
    
        sin = [] 
        for i in range(len(inflections_true) - 1):
    
            p1 = inflections_true[i]
            p2 = inflections_true[i+1]
            along_dist = vector.project(p2)-vector.project(p1)
            straight_dist = LineString([p1,p2]).length
            sin.append(along_dist / straight_dist)
        
        sin_reg = np.std(sin)
        sin_mean = np.mean(sin)
        # print(sin_mean, sin_reg)
    else:
        sinuosity = 0.00
        sin_reg   = 0.00

    # plt.show()
    if plot == True:
        plt.figure(figsize=  [6,4])
        plt.title(f'id: {ids}\nSinuosity: {np.round(sinuosity,8)}_{np.round(sin_reg,8)}\nWindow: {s}, width: {int(row.node_mwm)}, Degree: {degree}')
        plt.plot(*row.geometry.xy)
        plt.plot(*vector.xy)
        
        if len(inflections_true) > 1:
            plt.plot(*inf_line.xy)
        for p in inflections_true:
            plt.scatter(*p.xy)
        plt.axis('equal')
        plt.axis('off')
        plt.show()
    
    # return(sinuosity, sin_reg, new_df)
    return(sinuosity, sin_reg)


# ### Run the inflection point function

# In[381]:


x = [1,2]
y = [4,5]


np.cross(x, y)


# In[458]:


get_ipython().run_cell_magic('time', '', "# for i in [62263100071, 62263100081, 62263300011, 62263100051, 62261200961, 62263200191]:\n#     inflection_points([i], df_vector, True)\nD = df.copy()\nsins = []\nwidths = [5,10,15,20]\nwidths = [10, 15, 20]\ndegrees = [2, 4]\nids = df.iloc[0:5].reach_id.values\ndfs = []\nfor id in ids:\n    sins.append([])\n    for s in range(len(widths)):\n        for d in degrees:\n        \n            sin,reg,df_cross = inflection_points(id, D, True, widths[s],d , True)\n    \n            dfs.append(df_cross)\n            l = D[D.reach_id == id].iloc[0].geometry\n\n        # sins[i].append(sin)\n        # df_vector.loc[reach.index, 'sinuosity'] = sin\n        # df_vector.loc[reach.index, 'sin_reg'] = reg\n")


# In[433]:


tm = dfs[2].iloc[0]['target_min']
t  = dfs[2].iloc[0]['target']
tp = dfs[2].iloc[0]['target_plus']

print(t[0] - tm[0])
print(t[1] - tm[1])


# In[422]:


import matplotlib.colors as mcolors
colors = list(mcolors.CSS4_COLORS.keys())


# In[444]:


dfs[2]


# In[443]:


for i in range(4):
    A = dfs[2].iloc[i]
    
    print(np.cross(A.v1, A.v2), A.Cross)


# In[435]:


# for i in range(dfs[2].shape[0]):

for i in range(2):
    plt.scatter(*dfs[2].iloc[i].target, color = colors[i])
    plt.scatter(*dfs[2].iloc[i].target_min, color = colors[i])
    plt.scatter(*dfs[2].iloc[i].target_plus, color = colors[i])
plt.show()


# # __Clustering__

# In[433]:


from sklearn.cluster import KMeans
X = DF[['width_ratio_x', 'n_chan_max', 'sinuosity']].values
kmeans = KMeans(n_clusters=4, random_state=0, n_init="auto").fit(X)
kmeans.labels_
DF.loc[:,'cluster'] = kmeans.labels_

from mpl_toolkits import mplot3d

get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np
import matplotlib.pyplot as plt
fig = plt.figure()
ax = plt.axes(projection='3d')

# Data for three-dimensional scattered points
zdata = DF['width_ratio_x']
xdata = DF['sinuosity']
ydata = DF['n_chan_max']
ax.scatter3D(xdata, ydata, zdata, c=DF['cluster']);
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');

# df_vector_multi = df_main[df_main.n_chan_mod > 1]
# print(df_vector_multi.shape[0])

# X = df_vector_multi[['width_ratio', 'm_chn', 'sinuosity']].values
# kmeans = KMeans(n_clusters=2, random_state=0, n_init="auto").fit(X)
# kmeans.labels_
# df_vector_multi.loc[:,'cluster'] = kmeans.labels_


# fig = plt.figure()
# ax = plt.axes(projection='3d')

# # Data for three-dimensional scattered points
# zdata = df_vector_multi['width_ratio']
# xdata = df_vector_multi['sinuosity']
# ydata = df_vector_multi['m_chn']
# ax.scatter3D(xdata, ydata, zdata, c=df_vector_multi['cluster']);
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z');
# print(df_vector_multi[df_vector_multi.cluster == 1])


# # __Valley channel__
# - Calculate sinuosity
# - Valley slope

# #### Change multi lines to lines or return line

# In[6]:


def multiline_check(line):
    if line.geom_type == 'MultiLineString':
        line = shapely.ops.linemerge(line)
    
    if line.geom_type == 'MultiLineString':
        l = np.zeros(len(line.geoms))
        for i in range(len(line.geoms)):
            l[i] = line.geoms[i].length
        line = line.geoms[l.argmax()]
    
    return line


# #### Smooth line section with savitsky-golay first order smoothing

# In[7]:


def SG_smoothing(line, w, full = False):
    # print('SG_smoothing1', len(line.coords))
    if full == False:
        line = line.simplify(0.5, preserve_topology=True)
        # print(len(line.coords))
        segment_dis = 50
        max_segments = line.segmentize(1)
        line  = line.segmentize(segment_dis)
        ratio = int(len(max_segments.coords) / len(line.coords))
        
        w = int(w / ratio)
        # print(len(line.coords))
    else:
        line  = line.segmentize(2)

    # print('SG_smoothing2', len(line.coords))
    if w > len(line.coords):
        # print('SG_smoothing2.1', len(line.coords))
        w = len(line.coords)
    # print('SG_smoothing3', w, 1)
    XSG = sc.signal.savgol_filter(line.xy[0], w, 1)
    # print('SG_smoothing4')
    YSG = sc.signal.savgol_filter(line.xy[1], w, 1)
    # print('SG_smoothing5')
    
    SG  = LineString(list(zip(XSG, YSG)))
    # print('SG_smoothing6')
    return SG


# #### Find the offset width and smooth buffered line

# In[9]:


# one times
def find_offset(line, width, smoothing):
    # print('find_offset')
    
    buffer_init = int((width / 2))
    
    
    intersect_line_L = multiline_check(line.offset_curve(1  * buffer_init)) # L
    intersect_line_R = multiline_check(line.offset_curve(-1 * buffer_init)) # R
    L_length = intersect_line_L.length / line.length
    R_length = intersect_line_R.length / line.length
    poly = shapely.geometry.Polygon([*list(intersect_line_L.coords), 
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
        offset_L = shapely.geometry.LineString()
        offset_R = shapely.geometry.LineString()
    
    return offset_L, offset_R, buffer, L_length, R_length


# #### Calculate centerline from smoothed buffer (centerline_from_edges)

# In[11]:


# #OLD
# def centerline_from_edges(L, R, line, buffer, width, smoothing):
#     segments = np.linspace(start = 0,stop = 1, num = 300)
#     increase = 1.1
    
#     RL_check = False
    
#     c = 0
#     error = False
#     while RL_check == False:
#         # print(f'Counter RL_check: {c}')
#         Lof = L.parallel_offset(30)
#         Rof = R.parallel_offset(-30)

#         cent = []
#         RLs = []
       
#         # plt.figure()
#         # plt.title('centerline_from_edges')
#         # plt.plot(*L.xy)
#         # plt.plot(*R.xy)
#         # plt.plot(*line.xy)

#         if L.length > R.length:
#             p_line = L
#             d_line = R
            
#         else:
#             p_line = R
#             d_line = L
#         first_pd = d_line.interpolate(0, normalized = True)
#         first_pp = p_line.interpolate(0, normalized = True)        
#         RL = LineString([first_pd,first_pp])
#         RLs.append(RL)
#         cent.append(RL.interpolate(RL.length / 2))

        
#         for i in range(len(segments)):
    
#             # Lp = L.interpolate(segments[i], normalized = True)
#             # Rp = R.interpolate(segments[i], normalized = True)
            
#             # RL = LineString([Lp, Rp])
#             # RL_short = LineString(RL.segmentize(1).coords[1:-1])    
#             # RLs.append(RL)

            
#             Lp = p_line.interpolate(segments[i], normalized = True)
#             Rp = shapely.ops.nearest_points(Lp, d_line)[1]
#             RL = LineString([Lp, Rp])
#             RL_short = LineString(RL.segmentize(1).coords[1:-1])    
#             RLs.append(RL)
#             cent.append(RL.interpolate(RL.length / 2))
#             if RL.intersects(line):
                
#                 # cent.append(RL.interpolate(RL.length / 2))
#                 # plt.plot(*RL.xy, 'black', linewidth = 1)
#                 # print(RL_short.intersects(L), RL_short.intersects(R), i)
#                 # if (RL.intersects(Lof)):
#                 #     inter = RL.intersection(L)
#                 #     print(inter.geom_type)
#                 # cent.append(RL.interpolate(RL.length / 2))    
#                 if (RL_short.intersects(L)) | (RL_short.intersects(R)):
#                     # plt.plot(*RL.xy, 'black', linewidth = 1)
#                     RL_check = True
#                     # plt.plot(*RL_short.xy, 'red', linewidth = 1)
#                     # buffer *=increase
#                     # buff_line_L = multiline_check(line.parallel_offset(1   * buffer))
#                     # buff_line_R = multiline_check(line.parallel_offset(-1  * buffer))
                    
#                     # if (len(buff_line_L.coords) > 2) & (len(buff_line_R.coords) > 2):
#                     #     L = SG_smoothing(buff_line_L, width*smoothing)
#                     #     R = SG_smoothing(buff_line_R, width*smoothing)
#                     #     RL_check = False
#                     # else:
#                     #     # print('RL CHECK TRUE in RL int Lof or Rof')
#                     #     error    = True
#                     #     RL_check = True  

#                     # break;
#                     continue;
#                 else:
#                     cent.append(RL.interpolate(RL.length / 2))
#                     # print('RL CHECK TRUE')
#                     RL_check = True

#             #     if (RL.intersects(Lof)) | (RL.intersects(Rof)):
#             #         new_i = i
#             #         check = (RL.intersects(Lof)) | (RL.intersects(Rof))
#             #         while check:
#             #             new_i += 1
#             #             if new_i == len(segments):
#             #                 new_i = i
                            
                            
#             #             if RL.intersects(Lof):
#             #                 Lp = L.interpolate(segments[new_i], normalized = True)
#             #             else:
#             #                 Rp = R.interpolate(segments[new_i], normalized = True)
                        
#             #             RL = LineString([Lp, Rp])
#             #             check = (RL.intersects(Lof)) | (RL.intersects(Rof))
#             #             if ((new_i+1) == len(segments)) & (check == True):
#             #                 check = False
#             #         ####### ERROR POSSIBLE IF I TOO END INCREASE BUFFER SIZE!!!!!!!!!!
#             #  # print(i)
#             # plt.scatter(*Lp.xy)
#             # plt.scatter(*Rp.xy)
                
#         # plt.show()
#         ## UPDATE CENTERLINE
        
#         # RL_check = True
#         if c > 30:
#             print('Error in centerline from edges')
#             RL_check = True
#             error = True
        
#         c += 1
    
#     final_pd = d_line.interpolate(1, normalized = True)
#     final_pp = p_line.interpolate(1, normalized = True)
#     RL = LineString([final_pd,final_pp])
#     RLs.append(RL)
#     cent.append(RL.interpolate(RL.length / 2))
#     # print('centerline error: ', error)
#     if error == False:
#         # cent.append(final_p)
#         centLine = LineString(cent)

#         centLine = second_centerline(centLine, L, R, line)
#     else:
#         centLine = LineString()
#         RLs = []
#         L = LineString()
#         R = LineString()
#         buffer = np.nan

#     return centLine, RLs,L,R, buffer


# ##### NEW

# In[12]:


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

# In[18]:


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
# Function with different buffer sizes for each side

# In[13]:


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


# #### Clip raster image

# In[14]:


def raster_clip(dataset, clip):
    out_image, out_transform = mask(raster, clip, crop=True)
    out_meta = raster.meta
    
    out_meta.update({"driver": "GTiff",
          "height": out_image.shape[1],
          "width": out_image.shape[2],
          "transform": out_transform})
    return (out_image, out_meta)


# ## Valley Slope Download Images

# #### Get values along line from raster
# extract_slope_along_raster_line

# In[41]:


def extract_slope_along_raster_line(xarr, line, plot, ax, color, add_reg = True):
    profile = []
    dist    = [] 
    samples = 200

    for i in range(samples):
        point = line.interpolate(i  / samples - 1. , normalized=True)
        
        # access the nearest pixel in the xarray
        value = xarr.sel(x=point.x, y=point.y, method="nearest").data
        profile.append(value[0])
        dist.append(line.project(point))
    fit = np.polyfit(dist,profile,1)

    # print('min / max', np.array(profile).min(), np.array(profile).max(), np.min(profile) , np.max(profile))
    
    if plot == True:
        p = np.poly1d( fit)
        ax.plot(dist, profile, label = 'Profile', color = color)
        if add_reg == True:
            ax.plot(dist, p(dist), ls = '--', label = 'Regression', color = color)
        ax.set_xlabel('Distance (m)')
        ax.set_ylabel('Elevation (m)')
        ax.set_ylim([np.min(profile) * 0.95, np.max(profile) * 1.05])
        
    
    return fit[0], fit[1], dist, profile


# #### Retreive and save google earth image
# save_ee_image

# In[16]:


def save_ee_image(img, id, bounds, projection):
    fileName = directory + f'input_created/dem/{id}_dem.tif'

    if os.path.exists(fileName) == False:
        region = ee.Geometry.Rectangle([[bounds[0], bounds[1]], [bounds[2], bounds[3]]]);
    
        # Single-band GeoTIFF files wrapped in a zip file.
        url = img.getDownloadUrl({
            'name': f'{id}',
            'bands': ['dem'],
            'region': region,
            'format': "GEO_TIFF",
            'filePerBand':False
        })
        response = requests.get(url)
    
        
        fileName = directory + f'input_created/dem/{id}_dem.tif'
        with open(fileName, 'wb') as fd:
          fd.write(response.content)
    
    
    raster = rioxarray.open_rasterio(fileName)    
    raster = raster.rio.reproject(projection)

    return raster


# #### Calculate valley slope for center line and vallet line
# Valley_slope

# In[40]:

from rioxarray.merge import merge_arrays\ndef valley_slope_and_confinemend(DF, DF_node,projection, projection_ee, plot):\n    ee.Authenticate()\n    #Authenticate the Google earth engine with google account\n    ee.Initialize(project = \'ee-niekcde\')\n    img = ee.Image(\'MERIT/DEM/v1_0_3\')\n\n\n    df = DF.copy()\n    df_node = DF_node.copy()\n\n    \n    ids = df[df.vc.isna() == False].reach_id.values\n\n    width_multiplier = 5\n    # for id in ids:\n    for id in tqdm(ids,\n                 ncols = 50):\n        \n\n        row = df[df.reach_id == id]\n        node_row = df_node[df_node.reach_id == id]\n        \n    \n\n        if node_row.shape[0] == 0:\n            width = row.iloc[0].max_width\n        else:\n            width = node_row.max_width.mean()\n        \n        rowR = row.copy()\n        l = rowR.iloc[0].geometry\n        \n        # if l.buffer(width_multiplier * width).geom_type == \'MultiPolygon\':\n        #     l = SG_smoothing(l, 500, True).buffer(width_multiplier * width)\n        # else:\n        l = l.buffer(width_multiplier * width)\n       \n        rowR.loc[rowR.index, \'geometry\'] = l\n        rowR = rowR.to_crs(projection_ee)\n        \n        rowr = row.copy()\n        rowr = rowr.to_crs(projection_ee)\n        # rowV = row.set_geometry(\'vc\', crs = projection).to_crs(projection_ee)\n        \n        L = rowR.iloc[0].geometry.boundary\n        bounds = L.bounds\n\n        xmin = bounds[0]\n        xmax = bounds[2]\n        \n        ymin = bounds[1]\n        ymax = bounds[3]\n        \n        xdist = xmax - xmin\n        ydist = ymax - ymin\n\n        x1 = xmin\n        x2 = xmin + (xdist/2)\n        x3 = xmax\n        y1 = ymin\n        y2 = ymin + (ydist/2)\n        y3 = ymax\n        \n        \n        try:\n            print(\'1\')\n            raster = save_ee_image(img, id, bounds, projection)\n            min_val = raster[0,:,:].min()\n            max_val = raster[0,:,:].max()\n        except:\n            try:\n                print(\'2\')                \n                bounds1 = [x1, y1, x2, y3]\n                bounds2 = [x2, y1, x3, y3]\n                raster1 = save_ee_image(img, f\'{id}_1\', bounds1, projection)\n                raster2 = save_ee_image(img, f\'{id}_2\', bounds2, projection)\n                \n                max_val = np.array([raster1[0,:,:].max(), raster2[0,:,:].max()]).max()\n                min_val = np.array([raster1[0,:,:].min(), raster2[0,:,:].min()]).min()\n                # Merge/Mosaic multiple rasters using merge_arrays method of rioxarray\n                res = raster1.rio.resolution()\n                raster = merge_arrays(dataarrays = [raster1, raster2], res = (res[0], res[0]), crs=raster1.rio.crs,\n                                      vmin = min_val, vmax = max_val)\n                \n                f, [ax1,ax2,ax3] = plt.subplots(1,3, figsize = [12,7])\n                raster1[0,:,:].plot.imshow(ax = ax1, cmap = \'bwr\')\n                raster2[0,:,:].plot.imshow(ax = ax2, cmap = \'PRGn\')\n                raster[0,:,:].plot.imshow(ax = ax3, cmap = \'PRGn\')\n                ax3.plot(*row.iloc[0].geometry.xy)\n                plt.tight_layout()\n                plt.show()\n                \n                print(\'2.2\')\n            except:\n                print(\'except\')\n                try:\n                    print(\'4\')\n                    B1 = [x1, y1, x2, y2]\n                    B2 = [x2, y2, x3, y3]\n                    B3 = [x1, y2, x2, y3]\n                    B4 = [x2, y1, x3, y2]\n                    print(\'4.1\')\n                    raster1 = save_ee_image(img, f\'{id}_1\', B1, projection)\n                    raster2 = save_ee_image(img, f\'{id}_2\', B2, projection)\n                    raster3 = save_ee_image(img, f\'{id}_3\', B3, projection)\n                    raster4 = save_ee_image(img, f\'{id}_4\', B4, projection)\n                    max_val = np.array([raster1[0,:,:].max(), raster2[0,:,:].max(),\n                                        raster3[0,:,:].max(), raster4[0,:,:].max()]).max()\n                    min_val = np.array([raster1[0,:,:].min(), raster2[0,:,:].min(),\n                                        raster3[0,:,:].min(), raster4[0,:,:].min()]).min()\n                    print(\'4.2\')\n                    # Merge/Mosaic multiple rasters using merge_arrays method of rioxarray\n                    res1 = raster1.rio.resolution()\n                    \n                    raster = merge_arrays(dataarrays = [raster1, raster2, raster3, raster4],\n                                          res = (res1[0], res1[0]), crs=raster1.rio.crs)\n                    # Merge/Mosaic multiple rasters using merge_arrays method of rioxarray\n                    res = raster1.rio.resolution()\n                    raster = merge_arrays(dataarrays = [raster1, raster2, raster3, raster4],\n                                          res = (res[0], res[0]), crs=raster1.rio.crs)\n                    print(\'4.3\')\n                    print(raster[0,:,:].min(),raster[0,:,:].max())\n                    f, ax = plt.subplots(3,2, figsize = [12,7])\n                    raster1[0,:,:].plot.imshow(ax = ax[0,0], cmap = \'bwr\')\n                    raster2[0,:,:].plot.imshow(ax = ax[0,1], cmap = \'PRGn\')\n                    raster3[0,:,:].plot.imshow(ax = ax[1,0], cmap = \'PRGn\')\n                    raster4[0,:,:].plot.imshow(ax = ax[1,1], cmap = \'PRGn\')\n                    raster[0,:,:].plot.imshow(ax = ax[2,0], cmap = \'PRGn\', vmin = min_val, vmax = max_val)\n                    print(\'4.4\')\n                    # ax3.plot(*row.iloc[0].geometry.xy)\n                    plt.tight_layout()\n                    plt.show()\n                    print(\'4\')\n                except:\n                    print(\'still too big: \', id)\n        \n        if plot == True:\n            f, (ax1, ax2) = plt.subplots(1,2, figsize = (15,4))\n        \n            rs,_,_,_  = extract_slope_along_raster_line(raster, row.iloc[0].geometry, plot, ax1, \'red\') \n            if row.iloc[0].vc.is_empty:\n                vs = np.nan\n            else:\n                vs,_,_,_  = extract_slope_along_raster_line(raster, row.iloc[0].vc, plot, ax1, \'b\') \n        else:\n            rs,_,_,_  = extract_slope_along_raster_line(raster, row.iloc[0].geometry, plot, 1, \'red\') \n            if row.iloc[0].vc.is_empty:\n                vs = np.nan\n            else:\n                vs,_,_,_  = extract_slope_along_raster_line(raster, row.iloc[0].vc      , plot, 1, \'b\') \n\n        df.loc[row.index, \'RSlope\'] = rs\n        df.loc[row.index, \'VSlope\'] = vs\n\n    \n        if plot == True:\n            print(\'valley_slope3.6\')\n            # ax1.set_aspect(\'equal\', adjustable="datalim")\n            rs_string = f\'River Slope:  {np.round(rs, 6)}\'\n            vs_string = f\'Valley Slope: {np.round(vs, 6)}\'\n            \n            # print(\'valley_slope3.1\')\n            \n            ax1.set_title(f\'{rs_string.ljust(22)}\\n{vs_string.ljust(22)}\')\n            ax1.grid()\n            \n            raster[0,:,:].plot.imshow(ax = ax2, cmap = \'RdYlGn_r\',\n                                      vmin = min_val, vmax = max_val)\n            ax2.plot(*row.iloc[0].geometry.xy, color = \'r\', label = \'center\')\n            ax2.plot(*row.iloc[0].vc.xy, color = \'b\', label = \'Valley\')\n\n            ax2.plot(*l.boundary.xy, color = \'b\', label = \'Valley\')\n\n\n            print(raster.rio.bounds())\n            print(row.geometry.bounds)\n            # ax2.plot(*row.iloc[0].geometry.buffer(row.iloc[0].max_width).exterior.xy)\n            # ax2.plot(*row.iloc[0].vb_l.xy)\n            # ax2.plot(*row.iloc[0].vb_r.xy)\n\n            \n            ax2.set_title(f\'{id}\\nLengths: {int(row.iloc[0].geometry.length)}, {int(row.iloc[0].vc.length)}\')\n            ax2.legend()\n            # ax2.set_axis_off()\n            ax2.set_aspect(\'equal\', adjustable="datalim")\n            \n            plt.savefig(directory + f\'figures/valley_slope/{row.iloc[0].reach_id}.png\')\n            plt.show()\n        # Get confinemend measure\n        print(\'confinemend\')\n        width = node_row.max_width.mean()\n        Slope = get_orthogonals(row.iloc[0].geometry, width, row.iloc[0].reach_id, raster, projection, 25, False)\n        \n        df.loc[row.index, \'CSlopeL\']           = Slope[0]\n        df.loc[row.index, \'CSlopeR\']           = Slope[1]\n\n    return df\n    # return raster\n\n# W = valley_slope(G,df_node,projection, projection_ee, plot = False)\n')


# In[310]:


# CHECK Valley line
# max slope in plaats van mean slope
# kijk of je nan kan maken

31241200251
31241200171
31241300381


# In[341]:





# In[336]:


get_ipython().run_cell_magic('time', '', 'I = [31239000481\n,31238000331]\nI = [31239000481]\n\nD = df[df.reach_id.isin(I)]\nE = valley_slope_and_confinemend(D, df_node, projection, projection_ee, True)\nE\n')


# ## Valley slope convert shapely

# In[17]:


# def slopeFunction(line, projection, dem, plot = False, plotLabel = '', c = ''):
#     eeLine = ee.Geometry.LineString(list(line.coords) , proj = projection)

#     transect = geemap.extract_transect(
#         dem, eeLine, n_segments=10, reducer="mean", to_pandas=True)
#     transect = transect.dropna()
    
#     fit = np.polyfit(transect.iloc[:,0],transect.iloc[:,1],1)
#     p = np.poly1d( fit)

#     if plot == True:
#         plt.plot(transect.iloc[:,0],transect.iloc[:,1], label = plotLabel, color = c)
#         plt.plot(transect.iloc[:,0], p(transect.iloc[:,0]), ls = '--', label = f'{plotLabel} regres', color = c)
#         # plt.plot(transect.iloc[:,0],transect.iloc[:,1], transect.iloc[:,0], p(transect.iloc[:,0]), '-', label = f'{plotLabel} regres')
#     return fit[0]

# def valley_slope(DF, plot):
#     ee.Authenticate()
#     #Authenticate the Google earth engine with google account
#     ee.Initialize(project = 'ee-niekcde')
#     projection = 'EPSG:3857'

#     dem = ee.Image("MERIT/Hydro/v1_0_1").select('elv')
#     dem = dem.reproject(crs=projection, scale=100)


#     df = DF.copy()
#     ids = df[df.vc.isna() == False].reach_id.values
    
#     for i in ids:
#         # print(i)
#         r = df[df.reach_id == i]
#         # print(r)
#         row = r.iloc[0]
#         # print(row)
#         try:
#             SG = slopeFunction(row.geometry     , projection, dem, plot= plot, plotLabel = 'Centerline', c = 'red')
#             SC = slopeFunction(row.vc, projection, dem, plot= plot, plotLabel = 'Valley_line', c = 'blue')
#         except:
#             SG = np.nan
#             SC = np.nan
        
#         df.loc[r.index, 'RSlope']  = SG
#         df.loc[r.index, 'VSlope']  = SC
#         if (plot == True) & (np.isnan(SG) == False):
#             plt.title(f'{i}\nCenterline: {np.round(SG, 4)}\nValley Line: {np.round(SC, 6)}')
#             plt.legend()
#             plt.show()

#     return df


# ## Test valley center en slope

# #### 1

# In[172]:


ids = [df.iloc[3200 + 941 + 200].reach_id]
r = df[df.reach_id==ids[0]]
node_row = df_node[df_node.reach_id==ids[0]]
row = r.iloc[0]
line = row.geometry
list(line.coords)
buffer = 500

BLL = multiline_check(row.geometry.parallel_offset(1 * buffer))
BLR = multiline_check(row.geometry.parallel_offset(-1 * buffer))


plt.plot(*row.geometry.xy, linewidth = 5)

plt.plot(*BLL.xy, linewidth = 1)
plt.plot(*BLR.xy, linewidth = 1)

# print(row.geometry.length)


# #### 3

# In[1162]:


# test too short buffer zones
ids = D.reach_id.values
ids


# #### 3.1

# In[47]:


ids = [
31241200171]


# In[201]:


get_ipython().run_cell_magic('time', '', '# First need to open the files\n\nD = df.copy()\n\n# ids = [61650800941]\n# print(ids)\nvalley_center = []\nplot = True\n# for id in ids:\n# for id in [61622200841\n# ,61570100591\n# ,61550200211]:\nsinoss = []\n\ntmTotal = datetime.now()\nfor rep in range(1):\n    for id in ids[0:10]:\n        for s in [300]:\n        \n            tm = datetime.now()\n            print(id)\n            r = D[D.reach_id==id]\n            if s != 0:\n                r.loc[:,\'geometry\'] = SG_smoothing(r.iloc[0].geometry, s, True)\n    \n            row = r.iloc[0]\n            node_row = df_node[df_node.reach_id==id]\n        \n            \n            river_direction, off_L, off_R, RLs, buffer, L_length, R_length = river_dir(row,node_row, 50)   \n            # river_dir(row,node_row, 20)   \n            # river_direction = LineString()\n            if (river_direction.is_empty == False):\n                sinuosity = row.geometry.length / river_direction.length\n            if (river_direction.is_empty == False) & plot == True:\n            \n                bl = multiline_check(row.geometry.offset_curve(1  * int(node_row.max_width.mean()/2)))\n                br = multiline_check(row.geometry.offset_curve(-1 * int(node_row.max_width.mean()/2)))\n                bl_multi = row.geometry.offset_curve(1  * int(node_row.max_width.mean()/2))\n                br_multi = row.geometry.offset_curve(-1 * int(node_row.max_width.mean()/2))\n                \n                print(s, row.geometry.length, bl.length / row.geometry.length, br.length / row.geometry.length)\n                \n                bl_final = multiline_check(row.geometry.offset_curve(1 * buffer))\n                br_final = multiline_check(row.geometry.offset_curve(-1 * buffer))\n                \n                D.loc[r.index, \'vc\'] = river_direction\n                # D.loc[r.index, \'vb_l\'] = off_L\n                # D.loc[r.index, \'vb_r\'] = off_R\n                \n                sinuosity = row.geometry.length / river_direction.length\n                D.loc[r.index, \'vc_sinus\'] = sinuosity\n                \n                \n                f, ax = plt.subplots(figsize=(8, 6))\n                plt.axis(\'equal\')\n                plt.axis(\'off\')\n                if bl_multi.geom_type == \'MultiLineString\':\n                    for i in bl_multi.geoms:\n                        ax.plot(*i.xy,\'lawngreen\', linewidth = 4)\n                if br_multi.geom_type == \'MultiLineString\':\n                    for i in br_multi.geoms:\n                        ax.plot(*i.xy,\'lawngreen\', linewidth = 4)    \n                \n                ax.plot(*off_L.xy,\'r\', linewidth = 2, label = \'smoothed buffer\')\n                ax.plot(*off_R.xy,\'r\', linewidth = 2)\n            \n                ax.plot(*bl.xy,\'g\', linewidth = 2, label = \'river width\')\n                ax.plot(*br.xy,\'g\', linewidth = 2)\n            \n                ax.plot(*bl_final.xy,\'yellow\', linewidth = 2, label = \'unsmoothed buffer\')\n                ax.plot(*br_final.xy,\'yellow\', linewidth = 2)\n            \n                ax.plot(*row.geometry.xy, \'orange\', label = \'centerline\')\n\n                DFV = df_vector[df_vector.reach_id==id].iloc[0]\n                ax.plot(*DFV.geometry.xy, \'firebrick\', label = \'old centerline\')\n                \n                ax.plot(*river_direction.xy, \'b\', label = \'valley line\')\n        \n                ax.set_title(f"{id}\\nLength: {np.round(river_direction.length, 0)}, width: {np.round(row.node_mwm, 0)}\\nSinuosity: {np.round(sinuosity, 3)}")\n                \n                # for i in range(len(RLs)):\n                #     ax.plot(*RLs[i].xy, linewidth = 1, color = \'black\')\n                \n                plt.legend()\n                plt.savefig(directory + f\'figures/valley_center/{row.reach_id}.png\')\n                plt.show()\n            sinoss.append(sinuosity)\n            # print(datetime.now() - tm)\ndatetime.now() - tmTotal\n\n# 3.26 min # 0.85\n8.19 # 0.8\n# 3.3 min # 0.9\n')


# In[202]:


id = ids[0]
r = D[D.reach_id==id]
if s != 0:
    r.loc[:,'geometry'] = SG_smoothing(r.iloc[0].geometry, s, True)

row = r.iloc[0]
node_row = df_node[df_node.reach_id==id]


# In[175]:


int(node_row.max_width.mean()/2)


# In[181]:


r.rch_ids_up.values


# In[203]:


# bl_multi = row.geometry.offset_curve(1  * int(node_row.max_width.mean()/2))
# br_multi = row.geometry.offset_curve(-1 * int(node_row.max_width.mean()/2))

bl_multi = row.geometry.offset_curve(1  * 100)
br_multi = row.geometry.offset_curve(-1 * 100)

plt.plot(*row.geometry.xy)
# if bl_multi.geom_type == 'MultiLineString':
#     print('?')
#     for i in bl_multi.geoms:
#         plt.plot(*i.xy,'lawngreen', linewidth = 1)
# if br_multi.geom_type == 'MultiLineString':
#     for i in br_multi.geoms:
#         plt.plot(*i.xy,'lawngreen', linewidth = 1)   


# In[88]:


plt.plot(*Main.geometry.xy, 'r', linewidth = 5)
plt.plot(*Down.geometry.xy, linewidth = 3)
plt.plot(*orig.iloc[0].geometry.xy, 'g')


# In[58]:


a = 31241200171
b = 31241200591 
c = 31241200171 
d = 31241200591 
e = 31241200171

A == E


# #### slope calc test

# In[59]:


# ids_slopetest = df[df.vc_sinus.isna() == False].iloc[0:10].reach_id.values
# ids_slopetest = [44547900611
# ,44547900731
# ,44547900631]

# W = df[df.reach_id.isin(ids)]
W = D[D.reach_id.isin(ids)]
w = W.iloc[1]


# In[106]:


df.columns


# In[99]:


I = [31239000481
,31238000331]


# In[121]:


D = df[df.reach_id.isin(I)]
E = valley_slope_and_confinemend(D, df_node, projection, projection_ee, True)
E


# # __Confinemend__

# ## _Extrapolate line from point_
# def getExtrapoledLine

# In[19]:


def getExtrapoledLine(p1,p2, w, wm, single = False):
    'Creates a line extrapoled in p1->p2 direction'
    EXTRAPOL_RATIO = w*(wm/2) # 10 river width per side
    if single == True:
        a = p1
        b = (p1[0]+EXTRAPOL_RATIO*(p2[0]-p1[0]), p1[1]+EXTRAPOL_RATIO*(p2[1]-p1[1]) )
        # b = (p2[0]+EXTRAPOL_RATIO*(p1[0]-p2[0]), p2[1]+EXTRAPOL_RATIO*(p1[1]-p2[1]) )
    else:
        a = (p1[0]+EXTRAPOL_RATIO*(p2[0]-p1[0]), p1[1]+EXTRAPOL_RATIO*(p2[1]-p1[1]) )
        b = (p2[0]+EXTRAPOL_RATIO*(p1[0]-p2[0]), p2[1]+EXTRAPOL_RATIO*(p1[1]-p2[1]) )
    return LineString([a,b])


# ## calculate slope along orthogonal lines (double sided)
# def get_orthogonal

# In[20]:


def get_orthogonals(line,width, id,raster, projection, num_segments, plot = False):
    # start line after 300 meter
    # print('get_orthogonals')
    if line.length > 1e3:
        Skip = 200 / line.length
        
        segments = np.linspace(Skip,1 - Skip,num_segments)
    
        max_dif =[]
        
        LSlopes = []
        RSlopes = []
        # change to percentage interpolate
    
        # tm = datetime.now()
        # fileName = directory + f'input_created/dem/{id}_dem.tif'      
        # raster = rioxarray.open_rasterio(fileName) 
        # raster = raster.rio.reproject(projection)
        # print(datetime.now() - tm)
    
        # tm = datetime.now()
        for i in range(num_segments):
            min = plus =  0.0005
            if i == 0:
                min  = 0
            if i == (num_segments -1):
                plus = 0
            pp = line.interpolate(segments[i] + plus, normalized = True)
            pm = line.interpolate(segments[i] - min, normalized = True)
            l = LineString([pm, pp])
    
    
        
            ll = multiline_check(l.parallel_offset(1))
            lr = multiline_check(l.parallel_offset(-1))
    
            # plt.plot(*line.xy)
            # plt.plot(*l.xy)
            
            # plt.plot(*ll.xy)
            # plt.plot(*lr.xy)
            # plt.show()
            # tm = datetime.now()
            lol = getExtrapoledLine(list(l.interpolate(0.5, normalized=True).coords)[0],
                                   list(ll.interpolate(0.5, normalized=True).coords)[0], width*2, 2, single = True)
            lor = getExtrapoledLine(list(l.interpolate(0.5, normalized=True).coords)[0],
                                   list(lr.interpolate(0.5, normalized=True).coords)[0], width*2, 2, single = True)
    
            lol_min = LineString(lol.segmentize(1).coords[10::])
            lor_min = LineString(lor.segmentize(1).coords[10::])
            # print(datetime.now() - tm)
            # if (lol_min.intersects(line)) | (lor_min.intersects(line)):
            #     print(id)
            
            # tm = datetime.now()
            # p = shapely.ops.nearest_points(lol, line)[0]
            Sl, Il, Dl, Pl = extract_slope_along_raster_line(raster, lol, False, ax = '', color = 'red', add_reg = False) 
            Sr, Ir, Dr, Pr = extract_slope_along_raster_line(raster, lor, False, ax = '', color = 'red', add_reg = False) 
            # print(datetime.now() - tm)
            # if (plot == True) | (lol_min.intersects(line)) | (lor_min.intersects(line)):
            if (plot == True):
                f,(ax1,ax2) = plt.subplots(1,2, figsize=[10,5])
        
                P = Pl.copy()
                P.extend(Pr)
                max_dif.append(np.max(P) -np.min(P))
                ax1.plot(-1*np.array(Dl), Pl, 'cornflowerblue')
                ax1.plot(Dr, Pr, 'navy')
                
                ax1.vlines(0, 0,5000, color = 'b', linestyle = 'dashed')
                pol = np.poly1d([Sl, Il])
                ax1.plot(-1*np.array(Dl), pol(Dl), color = 'r', linestyle = (0, (3, 5, 1, 5)))
        
                pol = np.poly1d([Sr, Ir])
                ax1.plot(Dr, pol(Dr), color = 'r', linestyle = (0, (3, 5, 1, 5)))
                ax1.grid()
                ax1.set_ylim([np.min(P) * 0.95, np.max(P) * 1.05])
                ax1.set_xlabel('Length (m)')
                ax1.set_ylabel('height (m)')
                ax1.set_title(f'Profile with centerline position\nLeft Slope: {np.round(Sl, 4)}\nRight Slope: {np.round(Sr, 4)}')
    
    
                
                ax2.plot(*multiline_check(line.parallel_offset(width)).xy, 'lightcoral')
                ax2.plot(*multiline_check(line.parallel_offset(-width)).xy, 'lightcoral')
                ax2.set_aspect('equal')
                raster[0,:,:].plot.imshow(ax=ax2, cmap = 'RdYlGn_r')
                ax2.plot(*line.xy, 'r')
                ax2.plot(*lol.xy, 'cornflowerblue')
                ax2.plot(*lor.xy, 'navy')
                ax2.set_title(f'Elevation raster\nWidth: {int(width)}\nCross: {int(np.max(Dl) + np.max(Dl))}')
                
                
                f.suptitle(id)
                plt.tight_layout()
                plt.show()
    
            LSlopes.append(Sl)
            RSlopes.append(Sr)
        # print(datetime.now() - tm)
        # return (np.mean(LSlopes), np.std(LSlopes), np.mean(RSlopes), np.std(RSlopes))
        return (np.max(LSlopes), np.max(RSlopes))
    else:
        return (np.nan, np.nan)
# get_orthogonals_p(w.iloc[0].geometry, w.iloc[0].max_width, w.iloc[0].reach_id, projection)


# ## _Test confinmend_

# In[56]:


get_ipython().run_cell_magic('time', '', 'id = 44547900841\nid = 44306900411\nid = 44306600101\n\nids  = [44547900841, 44306900411, 44306600101]\nids = df.reach_id.values\n\nids = df[df.reach_id == df.iloc[1184].reach_id].reach_id.values\nids = df.iloc[1185::].reach_id.values\nMD = []\n\nT = 0\nfor id in tqdm(ids[0:3],\n                 ncols = 50):\n    # print(\'start: \', id)\n    w = df[df.reach_id == id]\n    wn = df_node[df_node.reach_id == id]\n\n    if w.iloc[0].geometry.length < 1e3:\n        print(\'Tooo short\')\n    else:\n        s = 200\n        w.loc[:,\'geometry\'] = SG_smoothing(w.iloc[0].geometry, s, True)\n            \n        \n        rowR = w.copy()\n        rowR.loc[rowR.index, \'geometry\'] = rowR.iloc[0].geometry.buffer(5 * rowR.iloc[0].max_width)\n        rowR = rowR.to_crs(projection_ee)\n        \n        # rowV = w.set_geometry(\'vc\', crs = projection).to_crs(projection_ee)\n        \n\n        L = rowR.iloc[0].geometry.boundary\n\n            \n            \n        bounds = L.bounds\n    \n        width = wn.max_width.mean()\n\n        # tm = float(datetime.now().strftime("%Y%m%d%H%M%S.%f"))\n        dem_slice = save_ee_image(img, w.iloc[0].reach_id, bounds, projection)\n        # A = (float(datetime.now().strftime("%Y%m%d%H%M%S.%f")) - tm)\n        # print(A)\n        # T += A\n        \n        for segms in [25]:\n            # print(\'segsm loop \', segms)\n        \n            Slope = get_orthogonals(w.iloc[0].geometry, width, w.iloc[0].reach_id, projection, segms, False)\n            MD.append(Slope)\n')


# # __Cut segments at manually added nodes__

# In[21]:


def remove_man_add_section(line, node, endP):
    P = line.interpolate(line.project(node))

    splits = shapely.ops.split(line, P.buffer(1))
    dist = np.zeros(len(splits.geoms))
    for s in range(len(splits.geoms)):
        dist[s] = endP.distance(splits.geoms[s])

    newLine = splits.geoms[dist.argmax()]

    return newLine
def remove_manual_add(DF, df_node):
    df = DF.copy()
    
    node_man = df_node[df_node.manual_add == 1].reach_id.unique()

    for id in node_man:
        # print(id)
        V_row = df[df.reach_id == id]
        # print(V_row)
        V = V_row.iloc[0]
        
        N = df_node[df_node.reach_id == id]
        NMA = N[N.manual_add == 1]
        NMN = N[N.manual_add == 0]
        
        N_S = N.node_id.min() in NMA.node_id.values
        N_E = N.node_id.max() in NMA.node_id.values
        NL = 'Not'
        if ((N_S == False) & (N_E  == False)) == False:
            if (N_S == True) & (N_E == True):
                if N.shape[0] == NMA.shape[0]:
                    df = df.drop(V_row.index, inplace = False)
                else:
                    node = N[N.node_id == NMN.node_id.min()].iloc[0].geometry
                    endP = Point(V.geometry.coords[0])
                    
                    NL = remove_man_add_section(V.geometry, node, endP)
    
                    node = N[N.node_id == NMN.node_id.max()].iloc[0].geometry
                    endP = Point(V.geometry.coords[-1])
                    
                    NL = remove_man_add_section(NL, node, endP)
            
            else:
                if (N_E == True):
    
                    node_id = NMN.node_id.max()
                    endP = Point(V.geometry.coords[-1])
                
                elif (N_S == True):
    
                    endP = Point(V.geometry.coords[0])
                    node_id = NMN.node_id.min()
    
            
                node = N[N.node_id == node_id].iloc[0].geometry
                NL = remove_man_add_section(V.geometry, node, endP)
            if NL != 'Not':
                df.loc[V_row.index, 'geometry'] = NL

    return df


# # __New segment definition__

# ## _MultiLine to LineString_

# In[22]:


def find_starting_coords(lineList):

    # Identify starting and ending points for two linestrings
    ds = []
    ind = []
    L1 = list(lineList[0].coords)
    L2 = list(lineList[1].coords)
    for i in [0,-1]:
        Pi = Point(L1[i])
        for j in [0,-1]:
            Pj = Point(L2[j])
            ds.append(Pi.distance(Pj))
            ind.append([i,j])
    indices = ind[np.array(ds).argmin()]

    # Combine coordinates for new linestring
    if  indices == [-1,-1]:
        L2.reverse()
        L1.extend(L2)
    elif indices == [0 ,-1]:
        L1.reverse()
        L2.reverse()
        L1.extend(L2)
    elif indices == [-1, 0]:
        L1.extend(L2)
    elif indices == [0 , 0]:
        L1.reverse()
        L1.extend(L2)
        
    linest = LineString(L1)
    return linest


# In[23]:


def multi_to_line(lineList):
    if len(lineList) == 2:
        # print('test')
        outLine = find_starting_coords(lineList)
    elif len(lineList) == 3:
        # print('test1')
        intermediateLine = find_starting_coords(lineList)
        
        lineList = [intermediateLine, lineList[2]]
        outLine = find_starting_coords(lineList)
    # plt.plot(*intermediateLine.xy, linewidth = 5)
    # plt.plot(*outLine.xy)
    return outLine


# ## _Calculate if segment is too short_

# In[24]:


def too_short(DF, DF_node, MW_ratio = 12, number_wave = 4):
    """Calculates if a reach length is too short based on the meander wavelength vs width of the river"""
    df = DF.copy()
    df_node = DF_node.copy()
    df_node = df_node[df_node.manual_add == 0]
    
    present = 'node_mwm' in df.columns
    if present == False:
        
        NMW = df_node.groupby('reach_id', as_index= False)['max_width'].mean()
        NMW = NMW.rename(columns={"max_width": "node_mwm"})
        df = pd.merge(df, NMW, on ='reach_id',how='left')

    
    df.loc[:, 'min_len'] = df.node_mwm * MW_ratio * number_wave
    df.loc[:, 'short'] = df.geometry.length < df.min_len
    return df


# ## _Find if the first or last node connects to given dataframe of nodes_

# In[25]:


def find_connected_side(node1, node2):
    nodeS = node2.geometry.distance(node1.geometry.iloc[0]).min()
    nodeE = node2.geometry.distance(node1.geometry.iloc[-1]).min()
    if nodeS > nodeE:
        rowConnect = [-1, -2]
    else:
        rowConnect = [0,1]
    return rowConnect


# ## _Get up and downstream id for each segment_
# - Manual addition means no up or downstream segment --> mostly in tributaries
# - Angle between segments is used to identify the correct upstream segment
# 
# Function connection:
# - Create_ud_id --> identify_up_down
#     - Looks if there is a viable up or downstream segment and then calls "get_az" to calc angle between segments
# - get_az --> find_connected_side, Afterwards calls azimuth
#     - find_connected_side: checks if the first or last node is connected to connected segment
#     - azimuth: checks angle between two nodes/points

# #### Calculate direction between two points

# In[26]:


def azimuth(point1, point2):
    '''azimuth between 2 shapely points (interval 0 - 360)'''
    angle = np.arctan2(point2.x - point1.x, point2.y - point1.y)
    return np.degrees(angle) if angle >= 0 else np.degrees(angle) + 360


def azimuth_coords(x1,y1,x2,y2):
    '''azimuth between 2 shapely points (interval 0 - 360)'''
    angle = np.arctan2(x2 - x1, y2 - y1)
    return np.degrees(angle) if angle >= 0 else np.degrees(angle) + 360


# #### Calculate azimuth between two segments (up or down)

# In[27]:


def get_az(id, id_ud, df , df_node):
    row_ud = df[df.reach_id == id_ud]
    if row_ud.shape[0] > 0:
        row_ud = row_ud.iloc[0]

        
        nodes    = df_node[df_node.reach_id == id]
        nodes_ud = df_node[df_node.reach_id == row_ud.reach_id]
        
        if (nodes_ud.shape[0] < 2) | (nodes.shape[0] < 2):
            az = 1e10
        else: 
            nodesSide   = find_connected_side(nodes, nodes_ud)
            nodesSideUd = find_connected_side(nodes_ud, nodes)
            
            az_M  = azimuth(nodes.iloc[nodesSide[1]].geometry,nodes.iloc[nodesSide[0]].geometry)
            az_UD = azimuth(nodes_ud.iloc[nodesSideUd[0]].geometry,nodes_ud.iloc[nodesSideUd[1]].geometry)
            # print(az_M, az_UD)
            if (270 < az_M < 360) & (0 < az_UD < 90):
                az = (360 - az_M) + az_UD
            elif (270 < az_UD < 360) & (0 < az_M < 90):
                az = (360 - az_UD) + az_M
            else:
                az = np.abs(az_M - az_UD)
            
            # az = np.abs(az_M - az_UD)

            
    else:
        az = 1e10
    return (az)


# #### Check Manual add

# In[28]:


def check_manual_add(df_node, id1, id2):
    nodes    = df_node[df_node.reach_id == id1]
    nodesUd  = df_node[df_node.reach_id == id2]
    nodesSide   = find_connected_side(nodes, nodesUd)
    nodesSideUd = find_connected_side(nodesUd, nodes) 

    manAdd   = nodes.iloc[nodesSide[0]].manual_add
    manAddUd = nodesUd.iloc[nodesSideUd[0]].manual_add
    if (manAdd != 0) | (manAddUd != 0):
        Add = True
    else:
        Add = False
    return Add


# #### Check confluence

# In[29]:


def confluence(id, df, ud_side):
    """Identify whether the down or upstream side is a confluence"""
    
    row = df[df.reach_id == id].iloc[0]
    if ud_side == 'up':
        rch_id = row.rch_id_up
    else:
        rch_id = row.rch_id_dn
    conf_type = 'No'
    conf_ids = []
    if len(rch_id) == 0:
        # print(f'No {ud_side} segment')
        conf = False    
    elif len(rch_id) == 11:
        
        if ud_side == 'dn':
            row_dn = df[df.reach_id == int(rch_id)].iloc[0]
            if (row_dn.rch_id_up == None):
                conf = False
                # print(f'No downstream confluence')
            else:
                if (len(row_dn.rch_id_up) <= 11):
                    conf = False
                    # print(f'No downstream confluence')
                else:
                    # print(f'Downstream Confluence')
                    conf = True
                    conf_type = 'DN'
                    conf_ids = list(map(int, row_dn.rch_id_up.split(' ')))
                
        else:
            # print(f'No upstream confluence')
            conf = False
    else:
        conf = True
        conf_type = 'Simple'
        conf_ids = list(map(int, rch_id.split(' ')))
        # print(f'{ud_side} Confluence')

    conf_ids.append(row.reach_id)
    return conf, conf_type, conf_ids


# In[30]:


def connect_in_conf(id, conf_ids, df,df_node):

    comb = list(itertools.combinations(conf_ids,2))

    angles = np.zeros(len(comb))
    for i in range(len(comb)):
        angles[i] = get_az(comb[i][0], comb[i][1], df, df_node)

    connected = comb[angles.argmin()]
    connected_check = id in connected
    return connected_check


# #### Identify up or down segment

# In[31]:


def identify_up_down(row, df,df_whole,df_node, updown):
    """selected Row, vector dataframe, unfiltered vector dataframe, node dataframe, up or down indication """
    max_angle = 60
    if updown =='up':
        n_rch = row.n_rch_up
        rch_id = row.rch_id_up
    else:
        n_rch = row.n_rch_dn
        rch_id = row.rch_id_dn   
    
    # print(f'identify_up_down: {updown}', row.reach_id)
    az_ret =np.nan
    if rch_id != None:
        # print(f'identify_up_down1')
        if len(rch_id) > 11: # check if SWORD UD reach contains more than one UD reach
            # print('identify_up_down1.1')
            
            # convert string of multiple reach id's into seperate ints
            ids_up = list(map(int, rch_id.split(' ')))
            rows_up = df[df.reach_id.isin(ids_up)] # Look for the id's in the dataframe

            # Update ids: remove ids that are not found in dataframe
            ids_up = rows_up.reach_id.values
            ids = []
            angles = []
            for ud in range(len(ids_up)):
                add      = check_manual_add(df_node, row.reach_id, ids_up[ud])
                angle_ud = get_az(row.reach_id, ids_up[ud], df, df_node)
                if (add == False):
                    ids.append(ids_up[ud])
                    angles.append(angle_ud)
            
            
            if len(ids) == 0: # if no ud segments found return NaN
                # print('identify_up_down1.1.1')
                rch_id = np.nan
            elif len(ids) == 1: # if one reach, select this reach if closest node is not mannually added
                # print('identify_up_down1.1.2')
                conf, conf_type, conf_ids = confluence(row.reach_id, df_whole, updown)
                # print(conf_ids)
                rch_id = ids[0]
                az_ret = angles[0]
                if (conf == True):
                    if conf_type == 'DN':
                        ##### NO DOWNSTREAM??
                        connected = connect_in_conf(row.reach_id, conf_ids, df,df_node)
                        if connected == False:
                            rch_id = np.nan
                            az_ret = np.nan                        
                    
                
            else:
                # print('identify_up_down1.1.3')
                
                # print(f'check for loop: {ids}')
                # print('identify_up_down1.1.3.4')
                # az = np.zeros(len(ids))
                
                # for ud in range(len(ids)):
                #     az[ud] = get_az(row.reach_id, ids[ud], df, df_node)
                angles   = np.array(angles)
                az_index = angles.argmin() 
                
                rch_id = ids[az_index]
                az_ret = angles[az_index]

                # print(f'return angle between segments: {az_ret}, {angles}')
        else:
            # print('identify_up_down1.2')
            rch_id = df[df.reach_id == int(rch_id)]
            if rch_id.shape[0] == 0:
                rch_id = np.nan
            else:
                # CHECK PRESENTS OF UD SEGMENT
                # Check manual add
                # CHECK IF CONNECTION TO UD SEGMENT IS CONNECTED RIVER OR MAIN CHANNEL/TRIBUTARY
                rch_id = rch_id.iloc[0].reach_id
                
                add1 = check_manual_add(df_node, row.reach_id, rch_id)
                angle_ud1 = get_az(row.reach_id, rch_id, df, df_node)
                # print(add1, angle_ud1)
                conf, conf_type,conf_ids = confluence(row.reach_id, df_vector_whole, updown)
                # print(conf_ids)
                az_ret = angle_ud1
                if (add1 == True) :
                    rch_id = np.nan
                    az_ret = np.nan
                elif (conf == True):
                    if conf_type == 'DN':
                        connected = connect_in_conf(row.reach_id, conf_ids, df,df_node)
                        if connected == False:
                            rch_id = np.nan
                            az_ret = np.nan
                    elif conf_type == 'Simple':
                        print('\n\n\n\n\n\n\n\n ERROR: conf type should not be possible.\n\n\n\n\n\n\n\n\n')
                        
    else:
        # print('identify_up_down2')
        rch_id = np.nan
    # print(f'final rch_id {updown}:', rch_id)
    return rch_id, az_ret

# A = 62221500101
# r = df_vector[df_vector.reach_id == A].iloc[0]
# rch_id_dn = identify_up_down(r,df_vector, df_node, 'dn' )
# rch_id_dn

# df_vector1[df_vector1.reach_id == 62221500266]


# #### Create up down id Function

# In[131]:


def create_ud_id(DF,df_whole, DF_node):
    df = DF.copy()
    df_node = DF_node.copy()
    
    df['up_segment'] = True
    df['dn_segment'] = True
    for i in tqdm(range(df.shape[0])):
    # for i in tqdm(range(1,2)):
        # print(i)
        r = df.iloc[i]
        
        rch_id_up, azup = identify_up_down(r,df,df_whole, df_node, 'up' )
        rch_id_dn, azdn = identify_up_down(r,df,df_whole, df_node, 'dn' )
        # print(f'reach_id:{r.reach_id}\n\t\tUp: {rch_id_up}\n\t\tDn: {rch_id_dn}')
        ind = df[df.reach_id == r.reach_id].index
        df.loc[ind, 'single_up'] = rch_id_up
        df.loc[ind, 'single_dn'] = rch_id_dn
        df.loc[ind, 'az_up'] = azup
        df.loc[ind, 'az_dn'] = azdn
    
        if np.isnan(rch_id_up):
            df.loc[ind, 'up_segment'] = False
        if np.isnan(rch_id_dn):
            df.loc[ind, 'dn_segment'] = False
    
    df['single_up'] = df.single_up.astype('Int64')
    df['single_dn'] = df.single_dn.astype('Int64')
    df['order'] = 0

    return df


# In[94]:


# 31241200591 31241200171 31241200591 31241200171
# 31241200161
check_ids = [31241200161, 31241200591, 31241200171]

df[df.reach_id == ids[0]]
A = df_vector_whole[df_vector_whole.reach_id.isin(check_ids)]
A


# ## Create new segments
# Combine segments as long as upstream of downstreamsegments are available. For each iteration check if the segment is too short

# #### Get up or downriver id 
# The id can be a number of orders removed from river segment

# In[159]:


def get_ud_id(row, df,dir, order):
    r = row
    for o in range(order):
        
        if dir == 'up':
            ud = r.single_up
        else:
            ud = r.single_dn
        # print(o, ud, dir)
        if (ud is pd.NA) == False:
            ud_row = df[df.reach_id == ud].iloc[0]
            r = ud_row
            id = r.reach_id
        else:
            id = None
            break;

    return id


# #### Connect up and downstream segments to segment

# In[244]:


df[(df.short) & ((df.up_segment == True) | (df.dn_segment == False))].shape[0] 
# df[(df.short)].shape[0]


# In[253]:


def connect_ud_geometries(DF, df_orig, order):
    df = DF.copy()
    
    ids = df[(df.short) & ((df.up_segment) | (df.dn_segment))].reach_id.values
    # print('connect_ud_geometries1', ids)
    for s in tqdm(ids,
                 ncols = 50):
        rw = df[df.reach_id == s]
        r = rw.iloc[0]
            
        l = []
        # print('connect_ud_geometries2')
        for dir in ['up', 'dn']:
            
            id = get_ud_id(r,df,dir, order)
            
            # print('connect_ud_geometries3', id, order)
            
            if (id != None) & (id != r.reach_id):
                rup = df_orig[df_orig.reach_id == id].iloc[0]
                l.append(rup.geometry)

                df.loc[rw.index, 'LCheck'] += rup.geometry.length
                df.loc[rw.index, f'rch_ids_{dir}'] += f'{str(rup.reach_id)} '
            else:
                df.loc[rw.index, f'{dir}_segment'] = False

            # Always add row geometry after the check for upstream segment
            if dir == 'up':
                l.append(r.geometry)
    
        if len(l) > 1:
            new_geom = shapely.geometry.MultiLineString(l)
            
            if new_geom.geom_type == 'MultiLineString':
                   
                new_geom = multi_to_line(list(new_geom.geoms))
            df.loc[rw.index, 'geometry'] = new_geom
            df.loc[rw.index, 'order']    = order
            
    
    return df


# #### Create the new segments

# In[255]:


def create_new_segments(DF, df_orig, DF_node):
    df      = DF.copy()
    df_orig = DF.copy()
    df_node = DF_node.copy()

    df.loc[:, 'LCheck'] = df.geometry.length
    df.loc[:, 'rch_ids_up'] = ''
    df.loc[:, 'rch_ids_dn'] = ''
    
    df = too_short(df,df_node, 12, 4)
    
    short = df[df.short == True].shape[0]
    ud = df[(df.short == True) & (df.up_segment == False) & (df.dn_segment == False)].shape[0]
    order = 1

    while short != ud:
        df = connect_ud_geometries(df, df_orig, order)
        df = too_short(df,df_node, 12, 4)
    
        short = df[df.short == True].shape[0]
        ud    = df[(df.short == True) & (df.up_segment == False) & (df.dn_segment == False)].shape[0]
        # print(short, ud)
    
        if order == 35:
            short = 1
            ud    = 1
        
        order += 1
        

    df.loc[:, 'reach_len'] = df.geometry.length
    return df


# ## RUN create upstream id and create new segments

# In[158]:


plt.plot(*B.iloc[0].geometry.xy, linewidth = 5, label = B.iloc[0].reach_id)
plt.plot(*B.iloc[1].geometry.xy, linewidth = 3, label = B.iloc[1].reach_id)
plt.plot(*B.iloc[2].geometry.xy, label = B.iloc[2].reach_id)
plt.legend()


# In[248]:


a = df[df.short]
a.iloc[2]
a = df_vector[df_vector.reach_id.isin([31105000551,31105000561])]
a


# In[254]:


BB = create_ud_id(a,df_vector_whole,df_node)
B = create_new_segments(BB, BB, df_node)


# In[249]:


B.rch_ids_up


# In[200]:


df[df.reach_id == ids[0]]


# In[632]:


df_upstream = create_ud_id(df_vector,df_vector_whole,df_node)


# In[638]:


## CHANGE FOR LOOP TO EXCLUDE SHORT SEGMENTS THAT DO NOT HAVE UP OR DOWNSTREAM SEGMENTS???
# If code is not finished after 35 iteration the code is stopped
df_new_segments = create_new_segments(df_upstream, df_upstream, df_node)


# # __Combine everything__

# In[ ]:





# In[36]:


def remove_file(path):
    if os.path.exists(path):
        os.remove(path)
    else:
        print(f"The file does not exist: {path}")

def remove_dems(df):
    ids = df.reach_id.values
    for id in ids:
        files = glob.glob(f'{directory}input_created/dem/{id}_*')
        for f in files:
            remove_file(f)


# In[37]:


directory = '/Users/niekcollotdescury/Desktop/PhD/RiverWidth/'

vector_files = glob.glob(directory + 'input/SWOT_vector/*shp')
node_files   = glob.glob(directory + 'input/SWOT_nodes/*shp')

projection    = 'EPSG:3857'
projection_ee = 'EPSG:4326'


# In[1157]:


#### create new done sa-61/sa-62/na 82/as-31/af-12

# CALC done


# In[76]:


c = 'as'
i = 31


# In[ ]:


# %%time
openFiles  = True
create_new = True
calc = True

# as 31
# as 45
# sa 61 
# sa 62

# for c, i in [['as', 31],['as', 45],['sa', 61],['sa', 62]]:
# save buffer and initial width values to check what initial buffer value needs to be set
init_buffer = []
init_width  = []
# ['na', 82]
for c, i in [['as', 31], ['as',45], ['af', 12], ['na', 82]]:
# for c, i in [['as', 31]]:
    if openFiles == True:
        # print('1')
        vector_save_file = directory +f'results/new_segments/vector/{c}_{i}_new_segments_vector.shp'
        node_save_file   = directory +f'results/new_segments/node/{c}_{i}_new_segments_node.shp' 
        # print('2')
        if (len(glob.glob(vector_save_file)) == 0) | (len(glob.glob(node_save_file)) == 0) | (create_new == True):
            
            print(f'Start Loop: {c}_{i}')
            vector_file = glob.glob(directory + f'input/SWOT_vector/{c}*{i}*shp')[0]
            node_file   = glob.glob(directory + f'input/SWOT_nodes/{c}*{i}*shp')[0]
            
            df_vector       = gpd.read_file(vector_file)
            df_vector_whole = df_vector.to_crs(projection)
            df_vector       = df_vector_whole[(df_vector_whole['type'] == 1) & (df_vector_whole['lakeflag'] == 0)]
            
            df_node = gpd.read_file(node_file)
            df_node = df_node.to_crs(projection)
            df_node = df_node[df_node.reach_id.isin(df_vector.reach_id)]
        
            print(f'remove manually added')
            # remove parts of segments with manually added nodes
            df_vector = remove_manual_add(df_vector,df_node)
        
            print('Create new segments')
            # create new segments from upstream and downstream nodes
            df_upstream     = create_ud_id(df_vector,df_vector_whole,df_node)
            df_new_segments = create_new_segments(df_upstream, df_upstream, df_node)
        
            df = df_new_segments.copy()
    
            # width ratio
            df = width_ratio(df, df_node)
    
            # Save files
            df.to_file(vector_save_file      
                       , driver='ESRI Shapefile', crs = projection)
            df_node.to_file(node_save_file      
                       , driver='ESRI Shapefile', crs = projection)
        else:
            df       = gpd.read_file(vector_save_file)
            df_node  = gpd.read_file(node_save_file)
        
    print('\nCalc new variables: width ratio, inflection sinuosity, regularity, valley sinuosity & valley slope\n')
    if calc == True:
        
    
        for id in tqdm(df.reach_id.values,
                 ncols = 50):
            # print(id)
        # for id in df.reach_id.values:
            reach = df[df.reach_id == id]
            
            sin10_2,reg10_2 = inflection_points(id, df, False, 10, 2)
            sin10_4,reg10_4 = inflection_points(id, df, False, 10, 4)
            sin15_2,reg15_2 = inflection_points(id, df, False, 15, 2)
            sin15_4,reg15_4 = inflection_points(id, df, False, 15, 4)
            sin20_2,reg20_2 = inflection_points(id, df, False, 20, 2)
            sin20_4,reg20_4 = inflection_points(id, df, False, 20, 4)
            
            df.loc[reach.index, 'sinInf10_2'] = sin10_2
            df.loc[reach.index, 'sinReg10_4'] = reg10_4
            df.loc[reach.index, 'sinInf15_2'] = sin15_2
            df.loc[reach.index, 'sinReg15_4'] = reg15_4
            df.loc[reach.index, 'sinInf20_2'] = sin20_2
            df.loc[reach.index, 'sinReg20_4'] = reg20_4


            # Smooth line before calculating valley line --> 
            # shapely.geometry.LineString.parallel_offset gives multiple linesegments otherwise
            # print(id,2)
            s = 300
            reach.loc[:,'geometry'] = SG_smoothing(reach.iloc[0].geometry, s, True)
            row = reach.iloc[0]
            node_row = df_node[df_node.reach_id == row.reach_id]
            
            
            river_direction, off_L, off_R, RLs, buffer,L_length, R_length = river_dir(row,node_row, 50)   

            init_buffer.append(buffer)
            init_width.append( node_row.max_width.mean())

            
            df.loc[reach.index, 'vc']           = river_direction
            df.loc[reach.index, 'Ll_check']     = L_length
            df.loc[reach.index, 'Rl_check']     = R_length
            df.loc[reach.index, 'vc']           = river_direction
            
            if river_direction.is_empty:
                df.loc[reach.index, 'vc_sinus'] = np.nan
            else:
                df.loc[reach.index, 'vc_sinus'] = row.geometry.length / river_direction.length
            

            
            # print(id,5)
        print('\nStart Valley slope calc\n')
        df = valley_slope_and_confinemend(df, df_node, projection, projection_ee, False)
        # remove_dems(df)
        
    
        print('start save files')
        
        # save all geometries to csv 
        df.to_csv(directory +f'results/all/{c}_{i}.csv'      , index=False)
        
        #save file without valley_center
        df_nv = df.copy()
        df_nv = df_nv.drop('vc', axis = 1)
        df_nv.to_file(directory +f'results/centerline/{c}_{i}_centerline.shp'      , driver='ESRI Shapefile', crs = projection)
    
        


# In[252]:


df


# # Open csv

# In[ ]:


############ check column 22

def open_mg_file(path, projection):
    df = pd.read_csv(path)
    print(df.columns)
    # decode geometry columns as strings back into shapely objects
    # for c in ["geometry","vc"]:
    for c in ["geometry","vc"]:
        df[c] = df[c].apply(shapely.wkt.loads)

    # finally reconstruct geodataframe
    gdf = gpd.GeoDataFrame(df, crs = projection)
    return gdf


# In[346]:


c = 'sa'
i = 62
fp = directory +f'results/all/{c}_{i}.csv'
df = open_mg_file(fp, projection)
df_node = gpd.read_file(directory + f'results/new_segments/node/{c}_{i}_new_segments_node.shp')
df.head()


# # testing

# In[967]:


##### Upstream downstream problems!!!!!!!!
#Wrong (cannont be both): 62221500101
# Normal:62222300091, 62234300151
# Large island: 62225700021
# Upstream identified as downstream: 62232900421, 62234200121


# In[75]:


print('\nStart Valley slope calc\n')
df = valley_slope_and_confinemend(df, df_node, projection, projection_ee, False)
# gdf = gpd.GeoDataFrame(df, crs = projection)
    
# Slope = get_orthogonals(w.iloc[0].geometry, width, w.iloc[0].reach_id, projection, segms, False)

print('start save files')
#save file without centerline
# df_nc = df.copy()
# df_nc = df_nc.drop('geometry', axis = 1)
# df_nc = df_nc.set_geometry('vc')
# df_nc.to_file(directory +f'results/valley/{c}_{i}_valley_center.shp'      , driver='ESRI Shapefile', crs = projection)

# save all geometries to csv 
df.to_csv(directory +f'results/all/{c}_{i}.csv'      , index=False)

#save file without valley_center
df_nv = df.copy()
df_nv = df_nv.drop('vc', axis = 1)
df_nv.to_file(directory +f'results/centerline/{c}_{i}_centerline.shp'      , driver='ESRI Shapefile', crs = projection)


# In[ ]:


I = [31239000481
,31238000331]


# In[ ]:


# gdf = gpd.GeoDataFrame(df, crs = projection)
    
# Slope = get_orthogonals(w.iloc[0].geometry, width, w.iloc[0].reach_id, projection, segms, False)

# print('start save files')
# #save file without centerline
# # df_nc = df.copy()
# # df_nc = df_nc.drop('geometry', axis = 1)
# # df_nc = df_nc.set_geometry('vc')
# # df_nc.to_file(directory +f'results/valley/{c}_{i}_valley_center.shp'      , driver='ESRI Shapefile', crs = projection)

# # save all geometries to csv 
# df.to_csv(directory +f'results/all/{c}_{i}.csv'      , index=False)

# #save file without valley_center
# df_nv = df.copy()
# df_nv = df_nv.drop('vc', axis = 1)
# df_nv.to_file(directory +f'results/centerline/{c}_{i}_centerline.shp'      , driver='ESRI Shapefile', crs = projection)


# In[51]:


df_vector


# In[47]:


print(f'Start: {c}_{i}')
vector_file = glob.glob(directory + f'input/SWOT_vector/{c}*{i}*shp')[0]
node_file   = glob.glob(directory + f'input/SWOT_nodes/{c}*{i}*shp')[0]
print('1')
df_vector       = gpd.read_file(vector_file)
print('2')
df_vector_whole = df_vector.to_crs(projection)
print('3')
df_vector       = df_vector_whole[(df_vector_whole['type'] == 1) & (df_vector_whole['lakeflag'] == 0)]


# In[48]:


get_ipython().run_cell_magic('time', '', "print('1', node_file)\n")


# In[49]:


get_ipython().run_cell_magic('time', '', "print('1', node_file)\ndf_node = gpd.read_file(node_file)\n")


# In[50]:


get_ipython().run_cell_magic('time', '', "print('2')\ndf_node = df_node.to_crs(projection)\nprint('3')\ndf_node = df_node[df_node.reach_id.isin(df_vector.reach_id)]\n")


# In[974]:


get_ipython().run_cell_magic('time', '', 'A = rowR.iloc[0].geometry\nA.geoms[0].exterior\nA.geoms[1].exterior\n\nB = shapely.ops.unary_union([A.geoms[0].exterior, A.geoms[1].exterior])\nshapely.line_merge(B).geom_type\n\na = row.iloc[0].geometry\nc = SG_smoothing(a, 500, True).buffer(10 * rowR.iloc[0].max_width)\nc.geom_type\n')


# In[650]:


A = df.copy()
# A.loc[:, 'B'] = init_buffer
df_novc = A[A.vc_sinus.isna() == True]


# In[651]:


df_novc


# In[710]:


plt.plot(*off_R.xy)
plt.plot(*off_L.xy)

poly = shapely.geometry.Polygon([*list(off_L.coords), 
                                     *list(off_R.coords)[::-1]])
plt.plot(*poly.exterior.xy)
D = gpd.GeoDataFrame({'id': 1,'geometry':[poly]})
D.to_file('/Users/niekcollotdescury/Downloads/polygon-centerline-master'      , driver='ESRI Shapefile', crs = projection)

poly.envelope.exterior.xy[0]


# In[761]:


from shapely.ops import unary_union
def fixedInterpolation(line, minx, miny, dist):
    count   = dist
    newline = []

    plt.figure()
    startpoint = [line.xy[0][0] - minx, line.xy[1][0] - miny]
    endpoint = [line.xy[0][-1] - minx, line.xy[1][-1] - miny]
    newline.append(startpoint)

    plt.scatter(*Point(startpoint).xy)
    while count < line.length:
        point = line.interpolate(count)
        newline.append([point.x - minx, point.y - miny])
        plt.scatter(*Point([point.x - minx, point.y - miny]).xy)
        count += dist

    newline.append(endpoint)
    plt.scatter(*Point(startpoint).xy)
    # plt.plot(*line.xy)
    plt.show()
    return newline

def densifyBorder(polygon, minx, miny, dist):
    if len(polygon.interiors) == 0:
        exterIN = LineString(polygon.exterior)
        points = fixedInterpolation(exterIN, minx, miny, dist)
    
    else:
        exterIN = LineString(polygon.exterior)
        points = fixedInterpolation(exterIN, minx, miny, dist)
    
        for j in range(len(polygon.interiors)):
            interIN = LineString(polygon.interiors[j])
            points += fixedInterpolation(interIN, minx, miny, dist)
    return points


def CENTERLINE(poly, dist):
    minx = int(min(inputGEOM.envelope.exterior.xy[0]))
    miny = int(min(inputGEOM.envelope.exterior.xy[1]))
    
    border = np.array(densifyBorder(inputGEOM, minx, miny, dist))

    
    vor = Voronoi(border)

    vertex = vor.vertices
    print(vertex)
    for i in vertex:
        print(i)
        ppp = Point(i)
        plt.scatter(*ppp.xy)
    plt.show()

    
    lst_lines = []
    for j, ridge in enumerate(vor.ridge_vertices):
        print(j, ridge)
        if -1 not in ridge:
            line = LineString([
                (vertex[ridge[0]][0] + minx, vertex[ridge[0]][1] + miny),
                (vertex[ridge[1]][0] + minx, vertex[ridge[1]][1] + miny)])
    
            if line.within(inputGEOM) and len(line.coords[0]) > 1:
                lst_lines.append(line)
    
    return unary_union(lst_lines)

dist = 600
inputGEOM = poly
CENTERLINE(poly, dist)


# In[748]:


plt.plot(*poly.envelope.exterior.xy)
plt.plot(*poly.exterior.xy)
L = LineString(poly.exterior)
plt.scatter(*L.interpolate(0).xy)
plt.scatter(*L.interpolate(1, normalized = True).xy)


# In[733]:


dist = 600
inputGEOM = poly
CENTERLINE(poly, dist)


# In[705]:


from shp2centerline import Shp2centerline

f = '/Users/niekcollotdescury/Downloads/polygon-centerline-master/polygon-centerline-master.shp'


# In[708]:


for d in [800]:
    print(d)
    A = Shp2centerline(f, 't.shp', d)
    q = gpd.read_file(directory + 'python/t.shp')
    Q = q.iloc[0].geometry
    print(Q.geom_type)
    Q = shapely.line_merge(Q)
    print(Q.geom_type)
    if Q.geom_type == 'LineString':
        plt.plot(*Q.xy)
        print(len(list(Q.segmentize(1).coords)))
        QS = SG_smoothing(Q, 500, True)
        plt.plot(*QS.xy)
    else:
        for i in Q.geoms:
            plt.plot(*i.xy)
    
    plt.plot(*off_R.xy)
    plt.plot(*off_L.xy)
    
    plt.show()


# # Analysis

# In[81]:


f, ax = plt.subplots()
gdf.vc_sinuosity.plot(kind = 'hist', ax = ax, bins = 60)
ax.set_title('valley center point sinuosity')


# In[82]:


f, ax = plt.subplots()
gdf.sinuosity.plot(kind = 'hist', ax = ax, bins = 60)
ax.set_title('inflection point sinuosity')


# In[83]:


f, ax = plt.subplots()
gdf.sin_reg.plot(kind = 'hist', ax = ax, bins = 60)
ax.set_title('inflection point sinuosity regularity')


# In[92]:


f, ax = plt.subplots()
gdf[['width_ratio', 'n_chan_mod']].plot(kind = 'hist', ax = ax, bins = 60, alpha = 0.5)
ax.set_title('width_ratio')
plt.xlim([0,10])

