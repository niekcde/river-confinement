#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 10:36:54 2024

@author: niekcollotdescury
"""

import os
import ee
import requests
import rioxarray

def save_ee_image(img, id, bounds, projection, directory):
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