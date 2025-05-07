# import ee
# ee.Authenticate()
# ee.Initialize('ee-niekcde')
def percDif(A, B):
    return abs(A-B) / ((A+B)/2)

import glob
Tiles = glob.glob('/scratch/6256481/results/new_segments/vector/*.shp')

import geopandas as gpd
import pandas as pd

import matplotlib.pyplot as plt

# dfWR = gpd.GeoDataFrame(columns = ['reach_id','Cont','Region','RL', 'RLC',
#                                     'PercDiff', 'geometry'])

# for i in range(len(Tiles)):
#     cont = Tiles[i][-29:-27] 
#     cont_reg = Tiles[i][-26:-24]
#     # print(cont,cont_reg, )
#     gdf = gpd.read_file(Tiles[i])

#     pdiff = percDif(gdf.reach_len,gdf.LCheck) 
#     WrongReaches = gdf[pdiff > 0.1]

#     pdiff = percDif(WrongReaches.reach_len,WrongReaches.LCheck)     
#     print(cont,cont_reg, WrongReaches.shape[0])
    

#     dfAdd= gpd.GeoDataFrame({'reach_id': WrongReaches.reach_id,'Cont':cont,'Region':cont_reg,
#                              'RL':WrongReaches.reach_len, 'RLC':WrongReaches.LCheck,
#                               'PercDiff': pdiff, 'geometry':WrongReaches.geometry})
#     dfWR = gpd.GeoDataFrame(pd.concat([dfWR, dfAdd]))

#     print(dfWR.shape[0])
# dfWR = dfWR.set_geometry('geometry')
# dfWR = dfWR.set_crs('EPSG:3857')

# print(dfWR.head())
# for i in range(dfWR.shape[0]):
#     plt.plot(*dfWR.iloc[i].geometry.xy)
#     plt.show()

# dfWR.to_csv('/scratch/6256481/results/WrongLengthReaches.csv')