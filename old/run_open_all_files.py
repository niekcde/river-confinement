import geopandas as gpd
import pandas as pd
import ast
import shapely
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import LineString
import scipy as sc


directory = '/scratch/6256481/'

import sys
sys.path.insert(0, directory + f'python/py_code/')


from open_files import open_all_files_multiproces
import glob
import warnings
warnings.filterwarnings("ignore")


sourceFiles  = np.sort(glob.glob(directory + 'input/SWOT_vector/*.shp'))
nodeFiles    = np.sort(glob.glob(directory + 'results/new_segments/node/*.shp'))
vectorFiles  = np.sort(glob.glob(directory + 'results/new_segments/vector/*.shp'))
files5       = np.sort(glob.glob(directory + 'results/all/*_5_5_0.csv'))
files10      = np.sort(glob.glob(directory + 'results/all/*10_5_0.csv'))
files15      = np.sort(glob.glob(directory + 'results/all/*15_5_0.csv'))

projection       = 'EPSG:3857'

# D = open_mg_file(files[0], ['geometry'], projection)
# listColumns  = ['apex5_0', 'ang5_0', 'outSlope', 'innSlope', 'CDists']
listColumns  = ['apex5_0', 'ang5_0', 'outSlope', 'innSlope', 'outPR',	'innPR', 'CDists', 'bendWidths']
pointColumns = ['infP5_0', 'apexP5_0']
# nestedLists  = ['CDists']
nestedLists  = ['outPR', 'innPR', 'CDists']
combMean     = [['outSlope', 'innSlope']]
combMeanName = ['mInnOutSlope']

# ###### import GLIM
GL = gpd.read_file(directory + 'input/GLIM.gdb')
GL = GL.to_crs(projection)
GL['GLGeometry'] = GL.geometry


size = 60
ERT  = 3
##Multiprocessing method:
F = np.sort(glob.glob(directory + f'results/all/*{size}_5_0.csv'))
sizes = []
for i in range(len(F)):
    with open(F[i]) as fp:
        c= 0
        for (count, _) in enumerate(fp, 1):
            c +=1
    sizes.append(c)
files = F[np.argsort(sizes)[::-1]]



def open_all_multi(F):
    A = open_all_files_multiproces(True, F,[-16, -11], ['geometry'], listColumns, nestedLists, combMean, combMeanName, directory, projection, GL , size, ERT)

from multiprocessing import Pool
if __name__ == '__main__':
    with Pool(10) as p:
        p.imap(open_all_multi, files)
        p.close()
        # # wait for all issued task to complete
        p.join()


from datetime import datetime
print('Done: ', datetime.now())