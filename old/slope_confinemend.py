#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 09:58:07 2024

@author: niekcollotdescury
"""

import ee
import numpy as np
import matplotlib.pyplot as plt


from save_ee_image import save_ee_image
from extract_slope_along_raster_line import extract_slope_along_raster_line
from get_orthogonals import get_orthogonals
    
from rioxarray.merge import merge_arrays


def valley_slope_and_confinemend(DF, DF_node,projection, projection_ee, plot, directory):
    # print('valley_slope_and_confinemend1')
    ee.Authenticate()
    
    # print('valley_slope_and_confinemend2')
    #Authenticate the Google earth engine with google account
    ee.Initialize(project = 'ee-niekcde')
    
    # print('valley_slope_and_confinemend3')
    img = ee.Image('MERIT/DEM/v1_0_3')
    
    # print('valley_slope_and_confinemend4')

    df = DF.copy()
    df_node = DF_node.copy()

    
    ids = df[df.vc.isna() == False].reach_id.values

    width_multiplier = 5
    # for id in ids:
    # for id in tqdm(ids,
    #              ncols = 50):
    for id in ids:
        

        row = df[df.reach_id == id]
        node_row = df_node[df_node.reach_id == id]
        
    

        if node_row.shape[0] == 0:
            width = row.iloc[0].max_width
        else:
            width = node_row.max_width.mean()
        
        rowR = row.copy()
        l = rowR.iloc[0].geometry
        
        # if l.buffer(width_multiplier * width).geom_type == 'MultiPolygon':
        #     l = SG_smoothing(l, 500, True).buffer(width_multiplier * width)
        # else:
        l = l.buffer(width_multiplier * width)
       
        rowR.loc[rowR.index, 'geometry'] = l
        rowR = rowR.to_crs(projection_ee)
        
        rowr = row.copy()
        rowr = rowr.to_crs(projection_ee)
        # rowV = row.set_geometry('vc', crs = projection).to_crs(projection_ee)
        
        L = rowR.iloc[0].geometry.boundary
        bounds = L.bounds

        xmin = bounds[0]
        xmax = bounds[2]
        
        ymin = bounds[1]
        ymax = bounds[3]
        
        xdist = xmax - xmin
        ydist = ymax - ymin

        x1 = xmin
        x2 = xmin + (xdist/2)
        x3 = xmax
        y1 = ymin
        y2 = ymin + (ydist/2)
        y3 = ymax
        
        # # print('1')
        # raster = save_ee_image(img, id, bounds, projection, directory)
        # min_val = raster[0,:,:].min()
        # max_val = raster[0,:,:].max()
        raster_check = 0
        try:
            # print('1')
            raster = save_ee_image(img, id, bounds, projection, directory)
            min_val = raster[0,:,:].min()
            max_val = raster[0,:,:].max()
        except:
            try:
                # print('2')                
                bounds1 = [x1, y1, x2, y3]
                bounds2 = [x2, y1, x3, y3]
                raster1 = save_ee_image(img, f'{id}_1', bounds1, projection, directory)
                raster2 = save_ee_image(img, f'{id}_2', bounds2, projection, directory)
                
                max_val = np.array([raster1[0,:,:].max(), raster2[0,:,:].max()]).max()
                min_val = np.array([raster1[0,:,:].min(), raster2[0,:,:].min()]).min()
                # Merge/Mosaic multiple rasters using merge_arrays method of rioxarray
                res = raster1.rio.resolution()
                raster = merge_arrays(dataarrays = [raster1, raster2], res = (res[0], res[0]), crs=raster1.rio.crs,
                                      vmin = min_val, vmax = max_val)
                
                f, [ax1,ax2,ax3] = plt.subplots(1,3, figsize = [12,7])
                raster1[0,:,:].plot.imshow(ax = ax1, cmap = 'bwr')
                raster2[0,:,:].plot.imshow(ax = ax2, cmap = 'PRGn')
                raster[0,:,:].plot.imshow(ax = ax3, cmap = 'PRGn')
                ax3.plot(*row.iloc[0].geometry.xy)
                plt.tight_layout()
                plt.show()
                
                # print('2.2')
            except:
                # print('except')
                try:
                    # print('4')
                    B1 = [x1, y1, x2, y2]
                    B2 = [x2, y2, x3, y3]
                    B3 = [x1, y2, x2, y3]
                    B4 = [x2, y1, x3, y2]
                    # print('4.1')
                    raster1 = save_ee_image(img, f'{id}_1', B1, projection, directory)
                    raster2 = save_ee_image(img, f'{id}_2', B2, projection, directory)
                    raster3 = save_ee_image(img, f'{id}_3', B3, projection, directory)
                    raster4 = save_ee_image(img, f'{id}_4', B4, projection)
                    max_val = np.array([raster1[0,:,:].max(), raster2[0,:,:].max(),
                                        raster3[0,:,:].max(), raster4[0,:,:].max()]).max()
                    min_val = np.array([raster1[0,:,:].min(), raster2[0,:,:].min(),
                                        raster3[0,:,:].min(), raster4[0,:,:].min()]).min()
                    # print('4.2')
                    # Merge/Mosaic multiple rasters using merge_arrays method of rioxarray
                    res1 = raster1.rio.resolution()
                    
                    raster = merge_arrays(dataarrays = [raster1, raster2, raster3, raster4],
                                          res = (res1[0], res1[0]), crs=raster1.rio.crs)
                    # Merge/Mosaic multiple rasters using merge_arrays method of rioxarray
                    res = raster1.rio.resolution()
                    raster = merge_arrays(dataarrays = [raster1, raster2, raster3, raster4],
                                          res = (res[0], res[0]), crs=raster1.rio.crs)
                    # print('4.3')
                    # print(raster[0,:,:].min(),raster[0,:,:].max())
                    
                    f, ax = plt.subplots(3,2, figsize = [12,7])
                    raster1[0,:,:].plot.imshow(ax = ax[0,0], cmap = 'bwr')
                    raster2[0,:,:].plot.imshow(ax = ax[0,1], cmap = 'PRGn')
                    raster3[0,:,:].plot.imshow(ax = ax[1,0], cmap = 'PRGn')
                    raster4[0,:,:].plot.imshow(ax = ax[1,1], cmap = 'PRGn')
                    raster[0,:,:].plot.imshow(ax = ax[2,0], cmap = 'PRGn', vmin = min_val, vmax = max_val)
                    # print('4.4')
                    # ax3.plot(*row.iloc[0].geometry.xy)
                    plt.tight_layout()
                    plt.show()
                    # print('4')
                except:
                    raster_check = 1
                    print('still too big: ', id)
        
        if plot == True:
            f, (ax1, ax2) = plt.subplots(1,2, figsize = (15,4))
        
            rs,_,_,_  = extract_slope_along_raster_line(raster, row.iloc[0].geometry, plot, ax1, 'red') 
            if row.iloc[0].vc.is_empty:
                vs = np.nan
            else:
                vs,_,_,_  = extract_slope_along_raster_line(raster, row.iloc[0].vc, plot, ax1, 'b') 
        elif raster_check == 0:
            rs,_,_,_  = extract_slope_along_raster_line(raster, row.iloc[0].geometry, plot, 1, 'red') 
            if row.iloc[0].vc.is_empty:
                vs = np.nan
            else:
                vs,_,_,_  = extract_slope_along_raster_line(raster, row.iloc[0].vc      , plot, 1, 'b') 
        else:
            rs = np.nan
            vs = np.nan

        df.loc[row.index, 'RSlope'] = rs
        df.loc[row.index, 'VSlope'] = vs

    
        if plot == True:
            print('valley_slope3.6')
            # ax1.set_aspect('equal', adjustable="datalim")
            rs_string = f'River Slope:  {np.round(rs, 6)}'
            vs_string = f'Valley Slope: {np.round(vs, 6)}'
            
            # print('valley_slope3.1')
            
            ax1.set_title(f'{rs_string.ljust(22)}\n{vs_string.ljust(22)}')
            ax1.grid()
            
            raster[0,:,:].plot.imshow(ax = ax2, cmap = 'RdYlGn_r',
                                      vmin = min_val, vmax = max_val)
            ax2.plot(*row.iloc[0].geometry.xy, color = 'r', label = 'center')
            ax2.plot(*row.iloc[0].vc.xy, color = 'b', label = 'Valley')

            ax2.plot(*l.boundary.xy, color = 'b', label = 'Valley')


            # print(raster.rio.bounds())
            # print(row.geometry.bounds)
            # ax2.plot(*row.iloc[0].geometry.buffer(row.iloc[0].max_width).exterior.xy)
            # ax2.plot(*row.iloc[0].vb_l.xy)
            # ax2.plot(*row.iloc[0].vb_r.xy)

            
            ax2.set_title(f'{id}\nLengths: {int(row.iloc[0].geometry.length)}, {int(row.iloc[0].vc.length)}')
            ax2.legend()
            # ax2.set_axis_off()
            ax2.set_aspect('equal', adjustable="datalim")
            
            plt.savefig(directory + f'figures/valley_slope/{row.iloc[0].reach_id}.png')
            plt.show()
        # Get confinemend measure
        # print('confinemend')
        if raster_check == 1:
            Slope = [np.nan, np.nan, np.nan]
            Slope4 = Slope6 = Slope8 = Slope
            print(Slope4, Slope6, Slope8)
        else:
            width = node_row.max_width.mean()
            Slope4 = get_orthogonals(row.iloc[0].geometry, width, row.iloc[0].reach_id, raster, projection, 2, 4, False)
            Slope6 = get_orthogonals(row.iloc[0].geometry, width, row.iloc[0].reach_id, raster, projection, 2, 6, False)
            Slope8 = get_orthogonals(row.iloc[0].geometry, width, row.iloc[0].reach_id, raster, projection, 2, 8, False)
        
        df.loc[row.index, 'CSlopeL4']           = Slope4[0]
        df.loc[row.index, 'CSlopeR4']           = Slope4[1]
        df.loc[row.index, 'CSlope4']            = Slope4[2]
        

        df.loc[row.index, 'CSlopeL6']           = Slope6[0]
        df.loc[row.index, 'CSlopeR6']           = Slope6[1]
        df.loc[row.index, 'CSlope6']            = Slope6[2]
        
        
        df.loc[row.index, 'CSlopeL8']           = Slope8[0]
        df.loc[row.index, 'CSlopeR8']           = Slope8[1]
        df.loc[row.index, 'CSlope8']            = Slope8[2]

    return df