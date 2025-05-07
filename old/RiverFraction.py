#!/usr/bin/env python
# coding: utf-8

# In[2229]:


# Add length of the segment!
# check intersection with larger river and stop segmentation a distance "width" of larger river from intersection
# check difference polygon and width --> to big a difference cut it off
# Bron picture van riv mask 


# # Package Import

# In[1]:


# Import Package 
import numpy as np

import re

from shapely import geometry, ops
from shapely.geometry import LineString, Point, Polygon
import shapely
from shapely.ops import nearest_points


import Read_inputs as ReadInputs
# import river_utils as ru

import scipy as sc
# from scipy import misc
from scipy.ndimage import binary_dilation, generate_binary_structure

import matplotlib.pyplot as plt

import geopandas as gpd
import pandas as pd
from multiprocessing import Pool
import rasterio
from rasterio import features

from pylab import *
import importlib

from osgeo import gdal

from tqdm import tqdm

import os

from datetime import datetime
max_gdal_cache_gb=64
gdal.SetCacheMax(int(max_gdal_cache_gb * 1e9))


# In[2231]:


importlib.reload(ReadInputs)


# In[2232]:


# import math

import matplotlib.pyplot as plt

import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle


# # Open or create functions

# ## Check the amount of unique water bodies in the tile

# In[4]:


def find_waterbodies(raster, keep_large = False):


    s8 = sc.ndimage.generate_binary_structure(2,2)
    label, _ = sc.ndimage.label(raster,structure = s8, output=np.uint64)
    
    waterbodies = np.unique(label)
    bodies = list()
    # for i in tqdm (waterbodies, 
    #            desc="Identify waterbodies…", 
    #            ascii=False, ncols=75):
    for i in waterbodies:
        if i != 0:
            labelc = label.copy()
            labelc[labelc != i] = 0
            d = raster.copy()
            d[labelc == 0] = 0
            bodies.append(d)
    #remove small bodies --> could be edge bodies
    if keep_large == True:
        large_bodies = []
        # for i in tqdm (range(len(bodies)), 
        #        desc="Keep large waterbodies…", 
        #        ascii=False, ncols=75):
        for i in range(len(bodies)):
            surface_area = np.count_nonzero(bodies[i])
            if surface_area > 5000:
                large_bodies.append(bodies[i])

        bodies = large_bodies
    
    return bodies


# ## Develop river and chanel mask

# In[5]:


def river_island_mask(data_points, im_number, dataset):
    print('river_island_mask1')
    river_mask_file  = directory + f'input_created/riv_mask/riv_mask_{im_number}.tif'
    island_mask_file = directory + f'input_created/isl_mask/isl_mask_{im_number}.tif'
    print('river_island_mask2')
    if (os.path.isfile(river_mask_file) == False) | (os.path.isfile(island_mask_file) == False):
        print('river_island_mask2.1')
        bodies = find_waterbodies(data_points, False)
        riv_mask_tile    = np.zeros([data_points.shape[0], data_points.shape[1]])
        island_mask_tile = np.zeros([data_points.shape[0], data_points.shape[1]])
        print('river_island_mask2.2')
        s8 = sc.ndimage.generate_binary_structure(2,2)
    
        label, _ = sc.ndimage.label(data_points,structure = s8, output=np.uint64)
        # waterbodies = np.unique(label)
        
        # for i in tqdm (range(len(bodies)), 
        #        desc="Create River and Island mask…", 
        #        ascii=False, ncols=75):
        for i in range(len(bodies)):
            
            chan_mask = bodies[i]
            chan_mask[chan_mask > 1] = 1
    
            
            operator = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]) #this is a 4-connect operator that is used for dilations and erosions
    
            #remove sections not connected to largest river section
            label, num_features = sc.ndimage.label(chan_mask,  structure = s8, output=np.uint64)
            if num_features > 1:
                labelhist, edges = np.histogram(label, bins=np.arange(label.max() + 2), density=False) 
                # labelhist_counter = labelhist[1:]
        
                maxid = np.where(labelhist[1:] == labelhist[1:].max())[0] + 1
                
                label[label != maxid[0]] = 0
                label[label > 0] = 1
                chan_mask = label.copy()
            
            # flipim is a version of the binary mask in which water=0 and land=1.  
            # Used to label and remove islands in the river.
            flipim = abs(1 - chan_mask)
            labelflip, _ = sc.ndimage.label(flipim, output=np.uint64)
            
            label      = 0
            chan_mask  = 0
            p          = 1
            dcnt       = 1 # adjust to get more dilation of the river??
        
            height = labelflip.shape[0]
            width  = labelflip.shape[1]
    
            # detect islands --> islands on the edge are not classified as an island
            while p == 1:
             
                s1 = labelflip[1, 0:width]
                s1s = np.histogram(s1, bins=np.arange(labelflip.max() + 2), density=False)
                if s1s[0][0] > 0:
                    s1s[0][0] = 0
                    
                s2 = labelflip[height -2, 0:width]
                s2s = np.histogram(s2, bins=np.arange(labelflip.max() + 2), density=False)
                if s2s[0][0] > 0:
                    s2s[0][0] = 0
    
                s3 = labelflip[0:height, 1]
                s3s = np.histogram(s3, bins=np.arange(labelflip.max() + 2), density=False)
                if s3s[0][0] > 0:
                    s3s[0][0] = 0
    
                    
                s4 = labelflip[0:height, width -2]
                s4s = np.histogram(s4, bins=np.arange(labelflip.max() + 2), density=False)
                if s4s[0][0] > 0:
                    s4s[0][0] = 0
    
                #Find unique values in the edges
                s1s0 = np.where(s1s[0] != 0)[0]
                s2s0 = np.where(s2s[0] != 0)[0]
                s3s0 = np.where(s3s[0] != 0)[0]
                s4s0 = np.where(s4s[0] != 0)[0]
                
                sidespoly = np.unique(np.concatenate([s1s0,s2s0,s3s0,s4s0]))
                sizepoly = sidespoly.shape[0]
    
                # change river parts to high values
                for y in range(sizepoly):
                    if labelflip[labelflip == sidespoly[y]].size != 0 and sidespoly[y] != 0:
                        labelflip[labelflip == sidespoly[y]] = 65535
                        
                #create Island mask
                island_mask = labelflip.copy()
                island_mask[(island_mask > 0) & (island_mask != 65535)] = 1
                island_mask[island_mask == 65535] = 0
                
                if labelflip[labelflip > 0].size != 0:          
                    labelflip[(labelflip > 0) & (labelflip != 65535)] = 0 # change islands to river value
                    labelflip[labelflip == 65535] = 1
                    
                    
                    riv_maskfl = labelflip
                    riv_mask = abs(1-labelflip)
                    
                    riv_maskd = binary_dilation(riv_mask, structure=operator).astype(riv_mask.dtype)
            
                    flipim = abs(1- riv_maskd )
                    labelflip, _ = sc.ndimage.label(flipim, output=np.uint64)
                    dcnt += 1
                    
                else:
                    p = 0
                            
                if dcnt > 1:
                    p = 0  
                        
            # Rivmask is the final river mask to be used in calculating the river center line.
            riv_mask = abs(1 - riv_maskfl)
        
            # Outputs the final river mask to a file.    
            riv_mask_tile    += riv_mask
            island_mask_tile += island_mask
            
        # print('create riv mask')
        with rasterio.open(river_mask_file, 'w',
                        driver = 'GTIFF',
                        height = riv_mask.shape[0],
                        width  = riv_mask.shape[1],
                        count = 1,
                        dtype = data_points.dtype,
                        crs = dataset.crs,
                        transform = dataset.transform) as dst: dst.write(riv_mask_tile, indexes = 1)
        #print('create island mask')
        with rasterio.open(island_mask_file, 'w',
                        driver = 'GTIFF',
                        height = island_mask.shape[0],
                        width  = island_mask.shape[1],
                        count = 1,
                        dtype = data_points.dtype,
                        crs = dataset.crs,
                        transform = dataset.transform) as dst: dst.write(island_mask_tile, indexes = 1)

    # riv_mask_dataset    = rasterio.open(river_mask_file)
    # island_mask_dataset = rasterio.open(island_mask_file)
    
    # riv_mask    = riv_mask_dataset.read(1)
    # island_mask = island_mask_dataset.read(1)
    # return (riv_mask, island_mask)


# ## Create riv boundary raster

# In[6]:


def create_riv_border(riv_mask, im_number, dir, data_points,dataset):

    file_path = dir + f'input_created/riv_mask_border/riv_mask_border_{im_number}.tif'

    if os.path.isfile(file_path) == False:

        operator = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]])
        imgbnd = riv_mask.copy()
        imgbnd[(riv_mask == 0) & (sc.ndimage.convolve(riv_mask, operator) > 0)] = 3
        imgbnd[(imgbnd == 0) & (sc.ndimage.convolve(imgbnd, operator) > 0)] = 4
        # imgbnd[(imgbnd == 0) & (sc.ndimage.convolve(imgbnd, operator) > 0)] = 4
        imgbnd[imgbnd < 4] = 0
        imgbnd[imgbnd == 4] = 1

        # Save raster
        with rasterio.open(file_path, 'w',
                            driver = 'GTIFF',
                            height = riv_mask.shape[0],
                            width  = riv_mask.shape[1],
                            count = 1,
                            dtype = data_points.dtype,
                            crs = dataset.crs,
                            transform = dataset.transform) as dst: dst.write(imgbnd, indexes = 1)


# ## Open River centerline vector and node files files
# Transform crs to crs of water mask

# In[7]:


def open_GWRL_SWOT(dataset, dir, im_number):

    river_vec_file = dir +f'input/GRWL_vector/{im_number}.shp'

    df_match = pd.read_csv(directory + 'input_created/GRWL_SWOT_match.txt')
    df_match.loc[df_match.Continent == 'NA_c', 'Continent'] = 'NA'
    df_match = df_match[df_match.Tile == im_number]
    
    
    SWOT_whole_file = dir+ f'input_created/SWOT_vector/SWOT_vector_{im_number}_whole.shp'
    SWOT_file       = dir+ f'input_created/SWOT_vector/SWOT_vector_{im_number}.shp'
    SWOT_nodes_file = dir+ f'input_created/SWOT_nodes/SWOT_nodes_{im_number}.shp'
    GWRL_file       = dir+ f'input_created/GWRL_vector/GWRL_vector_{im_number}.shp'

    if os.path.isfile(SWOT_file) == False:
        bounds  = dataset.bounds
        dataset_bounds_polygon = shapely.geometry.box(*bounds)
        
         # # GRWL reaches
        gdf_GWRL = gpd.read_file(river_vec_file)
        gdf_GWRL = gdf_GWRL.to_crs(dataset.crs)
        
        #Sword reaches
        SWOT_vector_files = []
        for i in range(df_match.shape[0]):
            match = df_match.iloc[i]
            swot_river_vector_file = directory +f'input/SWOT_vector/{match.Continent.lower()}_sword_reaches_hb{match.cont_number}_v16.shp'
            SWOT_vector_files.append(gpd.read_file(swot_river_vector_file))
        
        gdf_SWOT_whole = gpd.GeoDataFrame(pd.concat(SWOT_vector_files, ignore_index=True) )
        gdf_SWOT_whole = gdf_SWOT_whole.to_crs(dataset.crs)
        gdf_SWOT_whole = gdf_SWOT_whole[dataset_bounds_polygon.contains(gdf_SWOT_whole.geometry)].reset_index()
        
        # # Remove all different types and lakes. Keep only type River 
        gdf_SWOT_filtered = gdf_SWOT_whole[(gdf_SWOT_whole['type'] == 1) & (gdf_SWOT_whole['lakeflag'] == 0)]
        
        # # Split the reach_id in different pfafstetter parts
        gdf_SWOT = gdf_SWOT_filtered.copy()
        gdf_SWOT['pfaf']      = gdf_SWOT.reach_id.astype('str').map(lambda x: x[0:6]).astype('int')
        gdf_SWOT['r_number']  = gdf_SWOT.reach_id.astype('str').map(lambda x: x[6:10]).astype('int')
        gdf_SWOT['pfaf2']     = gdf_SWOT.reach_id.astype('str').map(lambda x: x[0:2]).astype('int')
    
        # Open SWOT Node file
        SWOT_node_files = []
        for i in range(df_match.shape[0]):
            match = df_match.iloc[i]
            swot_river_node_file = directory +f'input/SWOT_nodes/{match.Continent.lower()}_sword_nodes_hb{match.cont_number}_v16.shp'
            SWOT_node_files.append(gpd.read_file(swot_river_node_file))
        
        gdf_node = gpd.GeoDataFrame(pd.concat(SWOT_node_files, ignore_index=True) )
        gdf_node = gdf_node.to_crs(dataset.crs)
        gdf_node = gdf_node[dataset_bounds_polygon.contains(gdf_node.geometry)].reset_index()

        # add GWRL id to SWOT reach data 
        Mlines = []
        # for segid in tqdm (gdf_GWRL.segmentID.unique(), 
        #                desc="Identify segments…", 
        #                ascii=False, ncols=75):
        for segid in gdf_GWRL.segmentID.unique():

            
            segment = gdf_GWRL[gdf_GWRL.segmentID == segid]
            segment_multi_line = geometry.MultiLineString(list(segment.geometry.values))
            segment_line = ops.linemerge(segment_multi_line)
            Mlines.append(segment_line)
        
        segids = []
        # for i in tqdm (range(gdf_SWOT.shape[0]), 
        #                desc="Match GWRL id's with SWOT reach…", 
        #                ascii=False, ncols=75):
        for i in range(gdf_SWOT.shape[0]):
        
            reach = gdf_SWOT.iloc[i].geometry
        
            dists = reach.distance(Mlines)
            idmin = np.argwhere(dists == np.min(dists))
            segids.append(gdf_GWRL.segmentID.unique()[idmin[0][0]])
            
        gdf_SWOT['GWRL_ID'] = segids

    
        gdf_SWOT_whole.to_file(SWOT_whole_file, driver='ESRI Shapefile', crs = dataset.crs)
        gdf_SWOT.to_file(SWOT_file            , driver='ESRI Shapefile', crs = dataset.crs)
        gdf_node.to_file(SWOT_nodes_file      , driver='ESRI Shapefile', crs = dataset.crs)
        gdf_GWRL.to_file(GWRL_file            , driver='ESRI Shapefile', crs = dataset.crs)
        
    else:
        gdf_SWOT_whole = gpd.read_file(SWOT_whole_file)
        gdf_SWOT       = gpd.read_file(SWOT_file)
        gdf_node       = gpd.read_file(SWOT_nodes_file)
        gdf_GWRL       = gpd.read_file(GWRL_file)
    

    return (gdf_GWRL, gdf_SWOT_whole,gdf_SWOT,gdf_node)


# ## Polygonize raster tile
# - Polygonize raster tile
# - Check if polygonized file exists or create new file

# In[8]:


def polygonize_raster_tile(File,dataset, save=False,name = None):

    mask = None
    with rasterio.Env():
        with rasterio.open(File) as src:
            image = src.read(1) # first band
            results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v) 
            in enumerate(
                features.shapes(image, mask=mask, transform=src.transform)))
    # The result is a generator of GeoJSON features

    # change polygon to dataframe
    geoms = list(results)
    df = gpd.GeoDataFrame.from_features(geoms)
    

    border = df[df.raster_val == 1]
    border = shapely.unary_union(border.geometry)

    
    if save == True:
        df_border = gpd.GeoDataFrame({'geometry':[border]})
        df_border.to_file(name, driver='ESRI Shapefile', crs = dataset.crs)
    
    return border


# In[9]:


def open_create_polygonized(file,dataset, im_number,dir):

    file_polygon = dir + f'input_created/{file}_polygon/{file}_polygon_{im_number}.shp'

    if os.path.isfile(file_polygon):
        polygon = gpd.read_file(file_polygon).iloc[0,1]
    else:
        file_mask = dir + f'input_created/{file}/{file}_{im_number}.tif'
        polygon        = polygonize_raster_tile(file_mask, True, file_polygon)
    return (polygon)


# # Calculate functions

# ### parallel lines function
# Orthogonal_function

# In[10]:


def orthogonal_function(linestring, distance):
    """From a input shapely linestring with given distance, return right and left orthogonal points 
    only functions for straight linear lines"""
   
    plinel = linestring.parallel_offset(distance=distance, side = 'left')
    pliner = linestring.parallel_offset(distance=distance, side = 'right')
    
    return (plinel, pliner)


# ### Create line at the top or bottom of polygon
# polygon function in shapely can create strange shapes. This line can cut the top and buttom of of the polygon

# In[11]:


def border_line_side(X, center_line, side, disp):
    """X = Left or right parallel offset line, side is zero or 1 for the beginning or end of the line,
    center_line = the river center_line segment, disp = is distance from the centerline"""
    if side == 'start':
        X = X.coords[0]
        center = center_line.coords[0]
    elif side == 'mid':
        X = X.centroid.coords[0]
        center = center_line.centroid.coords[0]
    elif side == 'end':
        X = X.coords[-1]
        center = center_line.coords[-1]

    line = LineString([X, center])

    
    Xxsign = sign(X[0] - center[0])
    Xysign = sign(X[1] - center[1])

    Xx = X[0] + (disp * Xxsign)
    Xy = X[1] + (disp * Xysign)


    x = np.array([Xx, Xy])

    u = np.array(center)
    v = np.array(X)

    n = v - u
    n /= np.linalg.norm(n, 2)

    P = u + n*np.dot(x - u, n)

    newline = LineString(list(line.coords) + [P])
    
    return newline


# ### Function Adjust Raster polygon
# Reduce size of the raster tile. Set a box around river segment and reduce tile with that box

# In[12]:


def adjust_raster_polygon(poly, reach, reach_buffer, save = False, save_name = None):
    # tm = datetime.now() 

    #resize to reach area
    X = reach.geometry.iloc[0].xy[0]
    Y = reach.geometry.iloc[0].xy[1]
    Xmin = np.min(X)
    Xmax = np.max(X)
    Ymin = np.min(Y)
    Ymax = np.max(Y)
    
    # print(datetime.now() - tm)
    df_box = shapely.geometry.box(Xmin - reach_buffer, Ymin - reach_buffer, Xmax+ reach_buffer, Ymax+ reach_buffer)

    cutout = df_box.intersection(poly)
    
    return cutout

# river_mask_polygon_cutout        = adjust_raster_polygon(river_mask_polygon,df, 3e4)


# ### Find Intersect

# #### Get reach

# In[13]:


def get_reach(df, se, d = 100, get_longest = False):
    if se == 'start':
        s = d
    else:
        s = -d

    interPoint = df.geometry.interpolate(s).iloc[0]
    reach_splits = shapely.ops.split(df.iloc[0].geometry, interPoint.buffer(2))
    reach_lengths = np.zeros(3)
    # if split returns more that 3 segments there is a self intersect
    # below code remove self intersecting lines
    if len(list(reach_splits.geoms)) > 3:
        no_self_inter = []
        for i in range(len(list(reach_splits.geoms))):

            if (reach_splits.geoms[i].coords[0] ==reach_splits.geoms[i].coords[-1]) == False:
                no_self_inter.append(reach_splits.geoms[i])
        new_line = shapely.ops.linemerge(shapely.ops.unary_union(no_self_inter))

        reach_splits = shapely.ops.split(new_line, interPoint.buffer(1))
    
    
    for i in range(len(list(reach_splits.geoms))):
            reach_lengths[i] = reach_splits.geoms[i].length
    idmax = reach_lengths.argmax()
    idmin = reach_lengths.argmin()

    if ((d*2) > df.geometry.iloc[0].length) | (get_longest == True):
        reach_segment = list(reach_splits.geoms)[idmax]
    else:
        reach_segment = list(reach_splits.geoms)[np.where(reach_lengths == np.delete(reach_lengths, [idmax, idmin]))[0][0]]
    return (reach_segment)


# #### Line center function

# In[14]:


def line_center(line):
    """line is geometry object"""
    line_cent = line.interpolate((line.length / 2))
    return line_cent


# #### Get cent_border_line

# In[41]:


def cent_border_line(reaches, df_node,reach,reach_segment, centroid,river_border,start_side,side, width ,
                     reset_count_input = 0, displacement = 400):
    if reset_count_input != 0:
        reach_segment = get_reach(reach, start_side, displacement * (reset_count_input))
        R,L = orthogonal_function(reach_segment, 20)
        if side == 'R':
            centroid = line_center(R)
        else:
            centroid = line_center(L)
        reset_count = reset_count_input
    else:
        reset_count = 1

    buf_size = width / 2
    while_check = False
    R_inter_check,L_inter_check = orthogonal_function(reach_segment, 1)
    if side == 'R':
        inter_check = L_inter_check
    else:
        inter_check = R_inter_check
    while while_check == False:
        centroid_buf = centroid.buffer(buf_size)
        inters = centroid_buf.intersection(river_border)
        # reach_buf_inter = centroid_buf.intersection(reach_segment)
        # print(inters.is_empty)
        # plt.plot(*reach_buf_inter.xy)
        # plt.plot(*centroid_buf.exterior.xy)
        # plt.scatter(*centroid.xy)
        if inters.is_empty == False:
            near_point      = nearest_points(centroid, inters)[1]
            near_point_line = LineString([centroid,near_point])
            near_line       = LineString([line_center(reach_segment),near_point])
            # plt.title(f'{start_side}_{side}_{buf_size}')
            # plt.plot(*near_point.xy)
            # plt.plot(*near_point_line.xy)
            # plt.plot(*near_line.xy)
            # plt.show()
            if near_point_line.intersects(reaches):
                # print('\n\n\nINTERSECTS REACHES')
                if inters.geom_type == 'MultiPolygon':
                    # print('MULTIPOLYGON???????')
                    
                    inters_list = list(inters.geoms)
                    no_inter_list = []
                    c = []
                    for i in range(len(inters_list)):
                        N = nearest_points(centroid, inters_list[i])[1]
                        l = LineString([centroid,N])
                        L = LineString([line_center(reach_segment),N])
                        # print(l.intersects(reaches), L.buffer(1).intersection(reaches))

                        if (l.intersects(reaches) == False) & (L.buffer(1).intersection(reaches).geom_type 
                                                               == 'LineString') & (L.intersects(inter_check) == False):
                            no_inter_list.append(inters_list[i].boundary)
                            c.append(i)
                            break;
                        else:
                            # fig, ax = plt.subplots()
                            # ax.plot(*inters.geoms[i].exterior.xy, color = 'red')
                            # print('segmentize1')
                            if inters_list[i].area > 5e4:
                                Nb = nearest_points(centroid, inters_list[i])[1]
                                lb = LineString([centroid,Nb])
                                Lb = LineString([line_center(reach_segment),Nb])
                                # print(Nb)
                                ints = inters_list[i].boundary
                                for t in range(10):
                                    # print(t)
                                    remove_nearest = gpd.GeoDataFrame({'geometry':[Nb.buffer(50)]}, crs = reach.crs)
                                    dft = gpd.GeoDataFrame({'geometry':[ints]}, crs= reach.crs)
                                    dft2=dft.overlay(remove_nearest, how = 'difference')
                                    # dft2.plot(ax = ax, color = colors[t])
                                    ints = dft2.geometry.iloc[0]
                                    Nb = nearest_points(centroid, ints)[1]
                                    # print(Nb)
                                    lb = LineString([centroid,Nb])
                                    Lb = LineString([line_center(reach_segment),Nb])
                                    # ax.plot(*Nb.xy)
                                    # ax.plot(*lb.xy)
                                    # ax.plot(*Lb.xy)
                                    if (lb.intersects(reaches) == False) & (Lb.buffer(1).intersection(reaches).geom_type 
                                                                       == 'LineString') & (Lb.intersects(inter_check) == False):
                                        no_inter_list.append(ints)
                                        c.append(i)
                                        # print(t)
                                        break;
                                
                                # dft.plot(ax=ax, column = ids)
                            # print('segmentize3')
                            
                            # ax.set_title(inters_list[i].area)
                            # ax.plot(*reach_buf_inter.xy)
                            # ax.plot(*N.xy)
                            # ax.plot(*l.xy)
                            # ax.plot(*L.xy)
                            # plt.show()
                    
                    if len(no_inter_list) > 0:
                        # print(no_inter_list)
                        inters = shapely.ops.unary_union(no_inter_list)
                        near_point = nearest_points(centroid, inters)[1]
                        near_point_line = LineString([centroid,near_point])
                        near_line  = LineString([line_center(reach_segment),near_point])
                        
                    else:
                        near_point_line = None
                else:
                    near_point_line = None
            if (near_point_line != None) & (near_line.buffer(1).intersection(reaches).geom_type 
                                            == 'LineString') & (near_line.intersects(inter_check) == False):
                while_check = True                   
                final_l = LineString([line_center(reach_segment),near_point])


        if (buf_size > (20*width)) & (while_check == False):
            # print('reset: ', reset_count, while_check, start_side, side)
            reset_count += 1

            reach_segment = get_reach(reach, start_side, displacement * (reset_count))

            
            idnode = reach_segment.distance(df_node.geometry).argmin()
            width = df_node[df_node.index == idnode].iloc[0].width

            R,L = orthogonal_function(reach_segment, 20)
            if side == 'R':
                centroid = line_center(R)
            elif side == 'L':
                centroid = line_center(L)
            buf_size = width/2

            
            if (displacement*reset_count) > (reach.iloc[0].geometry.length*0.8):
                # print('displacement to far')
                while_check = True
                final_l = None

        buf_size += width / 4
    return (final_l, reset_count)


# #### Main find intersect function

# In[73]:


def find_intersect_lines(df, river_border, df_node, reaches):
    centroid_displacement = 300
    border_lines = []
    for side in ['start', 'end']:
        reach_segment = get_reach(df, side)
        
        idnode = reach_segment.distance(df_node.geometry).argmin()
        width  = df_node[df_node.index == idnode].iloc[0].width
        
        R,L = orthogonal_function(reach_segment, 20)
        R_centroid = line_center(R)
        L_centroid = line_center(L)

        lines = []
        reset_counter = []
        for line_side in ['R','L']:
            
            if line_side == 'R':
                bl,rc  = cent_border_line(reaches, df_node,df,reach_segment, R_centroid,river_border,side,line_side, width ,
                                         reset_count_input = 0, displacement = centroid_displacement)
                lines.append(bl)
                reset_counter.append(rc)
                # print(f'in IF statement: {side}_{line_side}_{bl}')
                # plt.plot(*bl.xy)
            elif line_side == 'L':
                if rc > 1:
                    bl,rc  = cent_border_line(reaches, df_node,df,reach_segment, L_centroid,river_border,side,line_side, width ,
                                         reset_count_input = reset_counter[0], displacement = centroid_displacement)
                else:
                    bl,rc  = cent_border_line(reaches, df_node,df,reach_segment, L_centroid,river_border,side,line_side, width ,
                                         reset_count_input = 0, displacement = centroid_displacement)
                # print(f'in IF statement: {side}_{line_side}_{bl}')
                # plt.plot(*bl.xy)
                # plt.show()
                lines.append(bl)
                reset_counter.append(rc)
                if (lines[0] != None) & (lines[1]!= None):
                    while reset_counter[0] != reset_counter[1]:
                        # print('reset unequal', reset_counter)
                        if reset_counter[0] < reset_counter[1]:
                            # print('R<L')
                            lines[0],reset_counter[0]  = cent_border_line(reaches, df_node,df,reach_segment, R_centroid,river_border,side,
                                                                          'R', width ,
                                                                          reset_count_input = reset_counter[1], 
                                                                          displacement = centroid_displacement)
                        else:
                            # print('R>L')
                            lines[1],reset_counter[1]  = cent_border_line(reaches, df_node,df,reach_segment, L_centroid,river_border,side,
                                                                          'L', width ,
                                                                          reset_count_input = reset_counter[0], 
                                                                          displacement = centroid_displacement)
            # print('Lines',lines, line_side, side)
        if None in lines:
            lines = None
            border_lines.append(lines)
        else:
            multi_line  = geometry.MultiLineString(lines)
            merged_line = ops.linemerge(multi_line)

            border_lines.append(merged_line)
    if None in border_lines:
        border_lines = None
    return border_lines


# ### Select_segment

# #### Split function

# In[74]:


def split_function(multi_poly, split, geomtype, select = False):
    if geomtype == 'Multi':

        for i in range(len(multi_poly.geoms)):
            # if (i == 1) | (i == 2):
               # plt.plot(*multi_poly.geoms[i].exterior.xy, label = i)
            S = multi_poly.geoms[i].intersects(split)
            # print(i)
            if S == True:
                # print('split', i)
                poly_id = i  
        # plt.legend()
        # plt.show()
        
        if select == True:
            poly = multi_poly.geoms[poly_id]
        elif select == False:
            poly = shapely.ops.split(multi_poly.geoms[poly_id], split)
    else:
        poly = shapely.ops.split(multi_poly, split)
    return poly


# #### Select segment

# In[75]:


def select_segment(df, river_mask, river_border, df_node, reaches):
    id = df.iloc[0,3]
    
    # print('SelectSegment1')
    # splits = find_intersect_lines(df, river_border, df_node, reaches)
    # print('Splits = : ', splits, len(splits))
    # print('end of find intersect 2')
    try:
        # print('SelectSegment1.1')
        splits = find_intersect_lines(df, river_border, df_node, reaches)
        # print('SelectSegment1.2')
        if splits != None:
            if (splits[0].geom_type == 'MultiLineString') or (splits[1].geom_type == 'MultiLineString'):
                # print('SelectSegment1.3')
                splits = None
                split_creation = 'False'
        else:
            splits = None
        # print('SelectSegment1.4')
        split_creation = 'True'
    except:
        # print('split_creation == False')
        split_creation = 'False'
        splits = None

    # print('SelectSegment2')
    if splits is not None:
        try:
            split_ids = []
            # print('SelectSegment2.1')
            if river_mask.geom_type == 'MultiPolygon':
                geomtype = 'Multi'
                for split in splits:
                    for i in range(len(river_mask.geoms)):
                        # print(f'for loop i: {i}')
                        S = river_mask.geoms[i].intersects(split)
                        if S == True:
                            # print(f'True i: {i}')
                            split_ids.append(i)
            else:
                geomtype ='Single'
            # print('SelectSegment2.2')
            # print(len(set(split_ids)))
            if len(set(split_ids)) <= 1:
                # print('same geom split')
                # print('SelectSegment2.3')
                
                first_poly  = split_function(river_mask, splits[0], geomtype, False)
                # plt.plot(*first_poly.geoms[0].exterior.xy, label = 'geom0')
                # plt.plot(*first_poly.geoms[1].exterior.xy, label = 'geom1')
                # plt.plot(*splits[0].xy,label = 'split0')
                # plt.plot(*splits[1].xy,label = 'split1')
                # plt.savefig(directory +'testfig.png')
                # plt.title('first_poly')
                # plt.legend()
                # plt.show()
                # print('SelectSegment2.4')
                second_poly = split_function(first_poly, splits[1], 'Multi', False)
                # print(list(first_poly.geoms))
                # plt.plot(*second_poly.geoms[0].exterior.xy, label = 'geom0')
                # plt.plot(*second_poly.geoms[1].exterior.xy, label = 'geom1')
                # plt.plot(*splits[0].xy,label = 'split0')
                # plt.plot(*splits[1].xy,label = 'split1')
                # plt.savefig(directory +'testfig.png')
                # plt.title('second_poly')
                # plt.legend()
                # plt.show()
                # print(second_poly.intersects(splits[0].buffer(1)))
                # print('SelectSegment2.5')
                poly = split_function(second_poly, splits[0].buffer(1), 'Multi', True) 
                # plt.plot(*poly.exterior.xy, label = 'geom0')

                # plt.plot(*splits[0].xy,label = 'split0')
                # plt.plot(*splits[1].xy,label = 'split1')
                # plt.savefig(directory +'testfig.png')
                # plt.title('final poly')
                # plt.legend()
                # plt.show()
                # print('SelectSegment2.6')
                # g = gpd.GeoDataFrame({'geometry':[poly],
                #                       'id':['reach_segment']})
                # g.plot(column = 'id', legend=False)
                # plt.title(f'single' )
                # plt.show()
                # print('SelectSegment2.7')
            else:
                # print('different geom split')
                first_poly   = split_function(river_mask, splits[0], geomtype, False)
                # print(list(first_poly.geoms))
                # plt.plot(*first_poly.geoms[0].exterior.xy, label = 'geom0')
                # plt.plot(*first_poly.geoms[1].exterior.xy, label = 'geom1')
                # plt.plot(*splits[0].xy,label = 'split0')
                # plt.plot(*splits[1].xy,label = 'split1')
                # plt.savefig(directory +'testfig.png')
                # plt.title('first_poly')
                # plt.legend()
                # plt.show()
                
                second_poly  = split_function(river_mask, splits[1], geomtype, False)
                # print(list(second_poly.geoms))
                # plt.plot(*second_poly.geoms[0].exterior.xy, label = 'geom0')
                # plt.plot(*second_poly.geoms[1].exterior.xy, label = 'geom1')
                # plt.plot(*splits[0].xy,label = 'split0')
                # plt.plot(*splits[1].xy,label = 'split1')
                # plt.savefig(directory +'testfig.png')
                # plt.title('first_poly')
                # plt.legend()
                # plt.show()
                # print('different geom split')
                reach = df.geometry.iloc[0]
                # short_line_p1 = shapely.ops.snap(reach.interpolate(200), reach, 0.5)
                # short_line_p2 = shapely.ops.snap(reach.interpolate(-200), reach, 0.5)
                # short_line_p  = shapely.geometry.MultiPoint([short_line_p1, short_line_p2])
                # short_line = shapely.ops.split(reach, short_line_p)
                # plt.plot(*short_line.xy)
                # dist = 1e20
                # pids = [0,0]
                fp = [p.intersection(reach).length for p in first_poly.geoms]
                sp = [p.intersection(reach).length for p in second_poly.geoms]
                # dm = np.argmax(d)
                # print(d, dm)
                # for fp in range(len(first_poly.geoms)):
                #     print(first_poly.geoms[fp].intersects(df.geometry.iloc[0]))
                #     print(first_poly.geoms[fp].intersection(reach).length)
                # for sp in range(len(second_poly.geoms)):
                #         print(second_poly.geoms[fp].intersects(df.geometry.iloc[0]))
                #         print(second_poly.geoms[fp].intersection(reach).length)
                #         d = first_poly.geoms[fp].distance(second_poly.geoms[sp])
                #         if d <= dist:
                #             dist = d
                #             pids = [fp,sp]
    
                poly = shapely.geometry.MultiPolygon([first_poly.geoms[np.argmax(fp)],second_poly.geoms[np.argmax(sp)]])
                # print('different geom split')
            poly_creation = 'True'
            
        except:
            # print('enter except')
            poly_creation = 'False'
            poly = None
    else:
        poly_creation = 'False'
        poly = None
    # print(poly_creation, split_creation)
    return (poly, splits, poly_creation, split_creation)


# ### Find areas with parallel reaches
# detect multilines

# In[76]:


def detect_multilines(df, border,dataset, buf = 30,remove = 1e5):
    df_buff = df.geometry.buffer(buf)
    df_buff = shapely.ops.unary_union(df_buff.geometry)

    
    poly_list = []
    # area_list = []
    for i in list(df_buff.geoms):
        all_internal_geoms = [geom for geom in i.interiors]
        poly_list.extend([shapely.geometry.Polygon(geom) for geom in all_internal_geoms])
    
    poly_list_large = []
    for p in poly_list:
        if (p.area > remove) & (p.intersects(border) ==False):
            poly_list_large.append(p)

    df_poly = gpd.GeoDataFrame({'geometry':poly_list_large}, crs=f'{dataset.crs}')

    df_poly_buff = df_poly.buffer(buf+15)
    poly_ids = []
    for r in range(df_poly_buff.shape[0]):
        
        inter = df[df_poly_buff.geometry[r].intersects(df.geometry)]
        inter_len = df_poly_buff.geometry[r].intersection(inter.geometry).length

        poly_ids.extend(inter[inter_len > 100].reach_id.values)

    df_return = pd.DataFrame({'id':poly_ids})
    
    return df_return


# In[ ]:





# ### Adjust reach up and down values

# In[77]:


def adjust_rch_up_dn(df, gdf):
    n_rch_up = df.n_rch_up
    n_rch_dn = df.n_rch_dn
    
    if (n_rch_up > 0) & (df.rch_id_up != None):
        if len(df.rch_id_up)>0:
            rch_id_up = np.array(df.rch_id_up.split(' ')).astype('int')
            rows_up   = gdf[(gdf.reach_id.isin(rch_id_up)) & (gdf.poly != 'False') & (gdf.split != 'False') & (gdf.in_river != 'False')]
            n_rch_up  = rows_up.shape[0]
        else:
            n_rch_up = 0
            rch_id_up = np.nan
    else:
        rch_id_up = np.nan
    if (n_rch_dn > 0) & (df.rch_id_dn != None):
        if len(df.rch_id_dn)>0:
            rch_id_dn = np.array(df.rch_id_dn.split(' ')).astype('int')
            rows_dn = gdf[(gdf.reach_id.isin(rch_id_dn)) & (gdf.poly != 'False') & (gdf.split != 'False') & (gdf.in_river != 'False')]
            n_rch_dn = rows_dn.shape[0]
        else:
            n_rch_dn = 0
            rch_id_dn = np.nan
    else:
        rch_id_dn = np.nan
    return (n_rch_up, n_rch_dn, rch_id_up, rch_id_dn)


# ### remove overlap

# In[78]:


def remove_overlap(df, df_poly, df_splits, df_whole):
    # print('\n\nremove_overlap')
    ids = []
    cut_polygons = []
    for r_id in df.reach_id.values:
        try:
            # print(f'\n\n{r_id}')
            
            reach = df_poly[df_poly.reach_id == r_id].iloc[0]
            reach_line = df[df.reach_id == r_id].iloc[0]
            
            # df_poly_r         = df_poly[df_poly.reach_id != r_id]
            df_splits_r       = df_splits[df_splits.reach_id != r_id]
            # df_r              = df[df.reach_id != r_id]
            df_whole_r        = df_whole[df_whole.reach_id != r_id]
            
            line_inters  = df_whole_r[reach.geometry.intersects(df_whole_r.geometry)]
            split_inters = df_splits_r[reach.geometry.intersects(df_splits_r.geometry)]
            # print(line_inters)
            line_collection = []
            # print(r_id)
            if line_inters.shape[0] > 0:
                lines = shapely.geometry.MultiLineString(line_inters.geometry.values)
    
                for l in lines.geoms:
                    # print('intersectlines: ', l.geom_type)
                    line_collection.append(l) 
        
                for l in range(split_inters.shape[0]):
                    line_collection.append(split_inters.iloc[l].geometry.geoms[0])
                    line_collection.append(split_inters.iloc[l].geometry.geoms[1])
                    # print('split1: ', split_inters.iloc[l].geometry.geoms[0].geom_type)
                    # print('split2: ', split_inters.iloc[l].geometry.geoms[1].geom_type)
                
                line_collection.append(reach.geometry.boundary)
                # print('boundary: ', reach.geometry.boundary.geom_type)
                # for i in range(len(line_collection)):
                # #     print(line_collection[i].geom_type)
                #     if line_collection[i].geom_type == 'MultiLineString':
                #         for g in line_collection[i].geoms:
                #             plt.plot(*g.xy)
                #         plt.show()
                    
                merged_lines = shapely.ops.linemerge(line_collection)
                border_lines = shapely.ops.unary_union(merged_lines)
                decomposition = shapely.ops.polygonize(border_lines)
    
                c = 0
    
                # plt.plot(*reach_line.geometry.xy)
                for d in decomposition:
                    # plt.plot(*d.exterior.xy, label = c)
                    if d.intersects(reach_line.geometry):
                        # print(c, d.intersects(reach.geometry))
                        cut_polygons.append(d)
                        ids.append(r_id)
                    c +=1
                #     plt.legend()
                # plt.show()
                    
            df_new_poly = gpd.GeoDataFrame({'id':ids, 'geometry':cut_polygons})
            new_poly = df_new_poly.groupby('id')['geometry'].apply(lambda x: shapely.ops.unary_union(x))     
            # print('new poly:', new_poly)
            for id in new_poly.index:
                # print(df[df.reach_id == id].index)
                index = df[df.reach_id == id].index[0]
                df_poly.loc[index, 'geometry'] = new_poly[id]
        except: 
            pass

    return df_poly


# ### Calculate fraction and number of islands

# In[79]:


def add_fraction_isles(gdf, gdf_poly, river_mask_polygon, island_mask_polygon):
    for i in range (gdf.shape[0]):
        # try:
            r = gdf.reach_id.values[i]
            reach_line = gdf[gdf.reach_id == r]
            reach_poly = gdf_poly[gdf_poly.reach_id == r].iloc[0]
        
            river_mask_polygon_cutout        = adjust_raster_polygon(river_mask_polygon       ,reach_line, 3e4)
            island_mask_polygon_cutout       = adjust_raster_polygon(island_mask_polygon      ,reach_line, 3e4)
            if reach_poly.poly == 'True':
                poly = reach_poly.geometry
                if poly.geom_type == 'Polygon':
                    reach_riv = river_mask_polygon_cutout.intersection(poly) 
                    riv_area = reach_riv.area
        
                    reach_isl = island_mask_polygon_cutout.intersection(poly)
                    isl_area = reach_isl.area
        
                    if reach_isl.geom_type == 'MultiPolygon':
                        num_isl = 0
                        for isl in reach_isl.geoms:
                            if isl.area > 500:
                                num_isl += 1
                    elif isl_area != 0.0:
                        num_isl = 1
                    else:
                        num_isl = 0
                    fraction = isl_area / riv_area
                    
                    # gdf.loc[reach_line.index, 'num_isl']      = num_isl
                    # gdf_poly.loc[reach_line.index, 'num_isl'] = num_isl
                    # gdf.loc[reach_line.index, 'frac']      = isl_area / riv_area
                    # gdf_poly.loc[reach_line.index, 'frac'] = isl_area / riv_area
                elif poly.geom_type == 'MultiPolygon':
                    riv_area = 0
                    isl_area = 0
                    num_isl  = 0
                    for g in poly.geoms:
                        reach_riv = river_mask_polygon_cutout.intersection(g) 
                        riv_area += reach_riv.area
            
                        reach_isl = island_mask_polygon_cutout.intersection(g)
                        isl_area += reach_isl.area
        
                        if reach_isl.geom_type == 'MultiPolygon':
                            for isl in reach_isl.geoms:
                                if isl.area > 500:
                                    num_isl += 1
                        elif isl_area != 0.0:
                            num_isl +=1 
                    fraction = isl_area / riv_area
                else:
                    num_isl  = np.nan
                    fraction = np.nan
                    

                gdf.loc[reach_line.index, 'num_isl']      = num_isl
                gdf_poly.loc[reach_line.index, 'num_isl'] = num_isl
                gdf.loc[reach_line.index, 'frac']      = fraction
                gdf_poly.loc[reach_line.index, 'frac'] = fraction
        # except:
        #     pass
        #     print('frac mistake')
            
    return (gdf, gdf_poly)


# # Open or create files Main
# - river mask
# - island mask
# - swot river vector and node files
# - river mask border
# 
#%%
def open_files(im_number):
    ##print('foldering')
    water_mask_file   = directory + 'input/GRWL_mask/' + im_number + '.tif'
    
    #print('open tile')
    dataset     = rasterio.open(water_mask_file)
    data_points = dataset.read(1)
    
    #print('Open or create river and island mask')
    # Open or create river and island mask
    riv_mask, island_mask = river_island_mask(data_points, im_number, dataset)
    
    
    #print('If not exists Create river border tif file')
    # If not exists Create river border tif file
    create_riv_border(riv_mask, im_number,directory, data_points,dataset)
    
    # #print('Open or create polygonized file')
    # # Open or create polygonized file
    # river_mask_polygon        = open_create_polygonized('riv_mask',        dataset, im_number,directory)
    # river_border_mask_polygon = open_create_polygonized('riv_mask_border', dataset, im_number,directory)
    # island_mask_polygon       = open_create_polygonized('isl_mask',        dataset, im_number,directory)
    
    # #print('open or create buffered river mask. To reduce the number of geometryies and connect broken small rivers')
    # # open or create buffered river mask. To reduce the number of geometryies and connect broken small rivers
    # riv_mask_polygon_buff_file = directory + f'input_created/riv_mask_polygon/riv_mask_polygon_buff_{im_number}.shp'
    # if os.path.isfile(riv_mask_polygon_buff_file) == False:
    #     river_mask_polygon_buff = river_mask_polygon.buffer(1)
    #     river_mask_polygon_buff = shapely.ops.unary_union(river_mask_polygon_buff)
    #     df_river_mask_polygon_buff = gpd.GeoDataFrame({'geometry':[river_mask_polygon_buff]}, crs=f'{dataset.crs}')
    #     df_river_mask_polygon_buff.to_file(directory + f'input_created/riv_mask_polygon/riv_mask_polygon_buff_{im_number}.shp'
    #                                        , driver='ESRI Shapefile')
    # riv_mask_polygon_buff = gpd.read_file(riv_mask_polygon_buff_file).iloc[0,1]


    # #print('Open and transform GWRL and SWOT files or open from created files')
    # # Open and transform GWRL and SWOT files or open from created files
    # gdf_GWRL, gdf_SWOT_whole,gdf_SWOT,gdf_SWOT_node = open_GWRL_SWOT(dataset, directory, im_number)
    
    # #print('Detect multi line sections in the tile')
    # # Detect multi line sections in the tile
    # multi_line_sections = detect_multilines(gdf_SWOT_whole,river_border_mask_polygon,dataset,remove = 1e5)  
    # return (dataset, data_points, 
    #         gdf_GWRL, gdf_SWOT_whole,gdf_SWOT,gdf_SWOT_node,
    #         river_mask_polygon, river_border_mask_polygon,island_mask_polygon, riv_mask_polygon_buff,
    #         multi_line_sections)
# In[80]:
tm = datetime.now()

def create_masks(im_number):
    ##print('foldering')
    print(im_number)
    water_mask_file   = directory + 'input/GRWL_mask/' + im_number + '.tif'
    print('create_masks1')
    #print('open tile')
    dataset     = rasterio.open(water_mask_file)
    data_points = dataset.read(1)
    print(data_points, dataset, im_number)
    print(type(data_points), type(dataset))
    print('create_masks2')
    #print('Open or create river and island mask')
    # Open or create river and island mask
    river_island_mask(data_points, im_number, dataset)
    print('create_masks3')
    
    #print('If not exists Create river border tif file')
    # If not exists Create river border tif file
    # create_riv_border(riv_mask, im_number,directory, data_points,dataset)
    print('create_masks4')
    return 1
import glob
directory         = '/Users/niekcollotdescury/Desktop/PhD/RiverWidth/'
tiles = [x[-8:-4] for x in glob.glob(directory + 'input/GRWL_mask/NN*.tif')]
tiles.extend([x[-8:-4] for x in glob.glob(directory + 'input/GRWL_mask/NO*.tif')])
tiles.extend([x[-8:-4] for x in glob.glob(directory + 'input/GRWL_mask/NP*.tif')])
# tiles.extend([x[-8:-4] for x in glob.glob(directory + 'input/GRWL_mask/NK*.tif')])
# tiles.extend([x[-8:-4] for x in glob.glob(directory + 'input/GRWL_mask/NL*.tif')])
# tiles.extend([x[-8:-4] for x in glob.glob(directory + 'input/GRWL_mask/NM*.tif')])

# tiles = [x[-8:-4] for x in glob.glob(directory + 'input/GRWL_mask/NM*.tif')]
# print(tiles)
create_masks('NN13')
# if __name__ == '__main__':
#     with Pool(2) as p:
#         A = p.map(create_masks, tiles)
# print( datetime.now() - tm)
  
#%%
f = directory + 'input/GRWL_mask/NM15.tif'
f = directory + f'input_created/riv_mask/riv_mask_NM16.tif'
rasterio.open(f)


# In[81]:


## %%time
#main df
def main_code(gdf_SWOT, river_mask_polygon, river_mask_polygon_buff, river_border_mask_polygon, 
              island_mask_polygon,
              gdf_SWOT_node, gdf_SWOT_whole,
              im_number,multi_line_sections,
              save):
    
    df = gdf_SWOT.copy()

    df = df[df.reach_id.isin(multi_line_sections.id) == False]

    
    gdf = df.copy()
    gdf['in_river'] = 'Unknown'
    gdf['frac']     = np.nan
    gdf['pfaf']     = gdf.reach_id.astype('str').map(lambda x: x[0:6]).astype('int')
    gdf['r_number'] = gdf.reach_id.astype('str').map(lambda x: x[6:10]).astype('int')
    gdf['poly']     = 'Unknown'
    gdf['split']    = 'Unknown'
    gdf['num_isl']  = np.nan
    #poly df
    gdf_poly = df.copy()
    gdf_poly['geometry'] = shapely.geometry.Polygon()
    gdf_poly['in_river'] = 'Unknown'
    gdf_poly['frac']     = np.nan
    gdf_poly['poly']     = 'Unknown'
    gdf_poly['split']    = 'Unknown'
    gdf_poly['num_isl']  = np.nan
    
    #splits df
    gdf_splits = df.copy()
    gdf_splits['geometry'] =  shapely.geometry.MultiLineString()
    gdf_splits['split'] = 'Unknown'
    
    reaches = shapely.ops.unary_union(gdf.geometry)
    
    # centerlines = []
    # counter = 0
    # print(gdf.shape[0])
    # next_ids = []
    # start_reaches = gdf[gdf.n_rch_up == 0].reach_id


    # for i in tqdm (range (gdf.shape[0]), 
    #                desc="Loading…", 
    #                ascii=False, ncols=75):
    for i in range (gdf.shape[0]):
        
        r = gdf.reach_id.values[i]
        df_reach = gdf[gdf.reach_id == r]
        # print(i, r)
        river_mask_polygon_cutout        = adjust_raster_polygon(river_mask_polygon       ,df_reach, 3e4)
        river_mask_polygon_buff_cutout   = adjust_raster_polygon(river_mask_polygon_buff  ,df_reach, 3e4)
        
        # island_mask_polygon_cutout       = adjust_raster_polygon(island_mask_polygon      ,df_reach, 3e4)
        river_border_mask_polygon_cutout = adjust_raster_polygon(river_border_mask_polygon,df_reach, 3e4)
        
        #change so that every reach that exists at the end of the reach is allowed!!!!!!!!!
        split_reach = shapely.ops.split(df_reach.geometry.iloc[0],river_border_mask_polygon_cutout)
        # plt.plot(*df_reach.iloc[0].geometry.xy, label = 'orig')
        if len(list(split_reach.geoms)) > 1:
            # print('split_reaches')
            num_inter = 0
            inter_segment = []
            for g in split_reach.geoms:
    
                if (g.intersects(river_mask_polygon_buff_cutout)) == True:
                    num_inter +=1
                    inter_segment.append(g)
            if num_inter == 1:
                df_reach.loc[df_reach.index, 'geometry'] = inter_segment[0]
            if num_inter > 1:
                s_len = 0
                s_id = 0
                # lens = []
                for s in range(num_inter):
                    if inter_segment[s].length > s_len:
                        s_len = inter_segment[s].length
                        s_id = s
                df_reach.loc[df_reach.index, 'geometry'] = inter_segment[s_id]
                if df_reach.iloc[0].geometry.touches(river_border_mask_polygon_cutout) == True:
                    reach_border_inter = df_reach.iloc[0].geometry.intersection(river_border_mask_polygon_cutout)
                    # plt.plot(*df_reach.iloc[0].geometry.xy, label = 'second')
                    # plt.scatter(*reach_border_inter.xy,c = 'red')
                    # plt.scatter(*Point(df_reach.iloc[0].geometry.coords[0]).xy)
                    # print(df_reach.iloc[0].geometry.length)
                    if Point(df_reach.iloc[0].geometry.coords[0]).distance(reach_border_inter) == 0:
                        se = 'start'
                        df_reach.loc[df_reach.index, 'geometry'] = get_reach(df_reach, se, d = 200, get_longest = True)
                        # print('yesstart')
                    else:
                        se = 'end'
                        df_reach.loc[df_reach.index, 'geometry'] = get_reach(df_reach, se, d = 200, get_longest = True)
                        # print('yesend')
                    # plt.plot(*df_reach.iloc[0].geometry.xy, label ='final')
        # plt.legend()
        # plt.show()
        centerline_in_river_complete = river_border_mask_polygon_cutout.crosses(df_reach.geometry.iloc[0])
        centerline_in_river_touches  = river_mask_polygon_cutout.intersects(df_reach.geometry.iloc[0])
        # print(centerline_in_river_complete, centerline_in_river_touches)
        if (centerline_in_river_complete == False) & (centerline_in_river_touches == True):
            gdf.loc[df_reach.index, 'in_river']      = 'True'
            gdf_poly.loc[df_reach.index, 'in_river'] = 'True'
    
            

            poly, splits, poly_creation, split_creation = select_segment(df_reach, 
                                                                river_mask_polygon_buff_cutout, 
                                                                river_border_mask_polygon_cutout, 
                                                                gdf_SWOT_node, reaches)
            # if poly.geom_type =="MultiPolygon":
            #     for p in poly.geoms: 
            #         plt.plot(*p.exterior.xy)
            # else:
            #    plt.plot(*poly.exterior.xy)
            # plt.plot(*df_reach.geometry.iloc[0].xy)
            
            # plt.plot(*splits[0].xy)
            # plt.plot(*splits[1].xy)
            # plt.show()
    
            gdf.loc[df_reach.index, 'poly']         = poly_creation
            gdf.loc[df_reach.index, 'split']        = split_creation
            gdf_poly.loc[df_reach.index, 'poly']    = poly_creation
            gdf_poly.loc[df_reach.index, 'split']   = split_creation
            gdf_splits.loc[df_reach.index, 'split'] = split_creation
            
        else:
            poly = None
            splits = None
            gdf.loc[df_reach.index, 'in_river']      = 'False'
            gdf_poly.loc[df_reach.index, 'in_river'] = 'False'
            gdf.loc[df_reach.index, 'poly']  = 'False'
            gdf.loc[df_reach.index, 'split'] = 'False'
            gdf_poly.loc[df_reach.index, 'poly']  = 'False'
            gdf_poly.loc[df_reach.index, 'split'] = 'False'
            gdf_splits.loc[df_reach.index, 'split'] = 'False'

        if splits is not None:
            g = shapely.geometry.MultiLineString(splits)
            gdf_splits.loc[df_reach.index, 'geometry'] = g
        else:
            gdf_splits.loc[df_reach.index, 'geometry'] = shapely.geometry.MultiLineString()
        if poly is not None:
            
            gdf_poly.loc[df_reach.index, 'geometry'] = poly
        else:
            gdf_poly.loc[df_reach.index, 'geometry'] = shapely.geometry.Polygon()
    # print(gdf_poly)
    try:
        gdf_poly = remove_overlap(gdf,gdf_poly,gdf_splits, gdf_SWOT_whole)
    except:
        #print('error in remove overlap: ', im_number)
        pass
    # print(gdf_poly)
    # try:
    gdf, gdf_poly = add_fraction_isles(gdf, gdf_poly, river_mask_polygon, island_mask_polygon)
    # except:
    #     print('frac calc: ', im_number)
    #     pass
    # print(gdf_poly)
    if save == True:
        gdf_poly.to_file(directory + f'result/poly_{im_number}.shp', driver='ESRI Shapefile')
        gdf.to_file(directory + f'result/vector_{im_number}.shp', driver='ESRI Shapefile')
        gdf_splits.to_file(directory + f'result/splits_{im_number}.shp', driver='ESRI Shapefile')

    return (gdf, gdf_poly, gdf_splits)


# # Main Code

# In[82]:


directory         = '/Users/niekcollotdescury/Desktop/PhD/RiverWidth/'
im_number         = 'SM19'
# (dataset, data_point, 
#  gdf_GWRL, gdf_SWOT_whole,gdf_SWOT,gdf_SWOT_node,
#  river_mask_polygon, river_border_mask_polygon,island_mask_polygon, river_mask_polygon_buff,
#  multi_line_sections) = open_files(directory, im_number)


# gdf, gdf_poly, gdf_split = main_code(gdf_SWOT, river_mask_polygon, river_mask_polygon_buff, river_border_mask_polygon, 
#                                      island_mask_polygon,
#                                      gdf_SWOT_node, gdf_SWOT_whole,
#                                      im_number, multi_line_sections,
#                                      save = True)
directory         = '/Users/niekcollotdescury/Desktop/PhD/RiverWidth/'
tiles = ['SM19','SM18','NG43']
def run_code(im_number):
    (dataset, data_point, 
     gdf_GWRL, gdf_SWOT_whole,gdf_SWOT,gdf_SWOT_node,
     river_mask_polygon, river_border_mask_polygon,island_mask_polygon, river_mask_polygon_buff,
     multi_line_sections) = open_files(directory, im_number)


    gdf, gdf_poly, gdf_split = main_code(gdf_SWOT, river_mask_polygon, river_mask_polygon_buff, river_border_mask_polygon, 
                                     island_mask_polygon,
                                     gdf_SWOT_node, gdf_SWOT_whole,
                                     im_number, multi_line_sections,
                                     save = True)
    return (gdf, gdf_poly, gdf_split)

# gdf, gdf_poly, gdf_split = run_code('SM19')

# if __name__ == '__main__':
#     with Pool(3) as p:
#         A = p.map(run_code, tiles)


# In[83]:


# Finished runs:
 # - SM18
    # 66120000451 problem with the split function --> the splits cross creating a smaller section
 # - NG43 --> frac calc error
 # - NG44
 # - NG45 Error
 # - NG46 --> frac calc error
 # - SB18
 # - SB19
 # - SB20
 # - SB21
# directory         = '/Users/niekcollotdescury/Desktop/PhD/RiverWidth/'
# tiles = ['SM18','NG43','NG44','NG45','NG46','SB18','SB19','SB20','SB21']
# tiles = ['SM19','SM18','NG43','NG44','NG45','NG46','SB18','SB19','SB20','SB21']
# for t in tiles:
    
#     df = gpd.read_file(directory + f'result/vector_{t}.shp')
#     print(t, df.shape[0], df[df.frac.isna()].shape[0], df[(df.in_river == 'True') & (df.poly == 'True') & (df.split == 'True')].shape[0])

