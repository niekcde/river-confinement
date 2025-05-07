import geopandas as gpd
import pandas as pd

import matplotlib.pyplot as plt

from support import open_mg_file

import glob
directory = '/scratch/6256481/'

FC = glob.glob(directory + 'results/all/?????.csv')
FA = glob.glob(directory + 'results/all/*conf*')
FC.sort()
FA.sort()

fc = [f[-9:-4] for f in FC]


fa = [f[-14:-9] for f in FA]


failed = []
for f in range(len(fc)):

    if fc[f] not in fa:
        failed.append(f)

files = glob.glob(directory + 'input/SWOT_vector/*.shp')
files.sort()

c = [ FC[i] for i in failed]
# print(len(fa), fa)

# print(glob.glob(directory + 'results/all/*conf*'))
import numpy as np
print(fc.index('na_74'))


from datetime import datetime
c = 'B'
i = 1
code_failure  = ''

code_failure += 'inflection_point failure'
code_failure += '\nconfiment failure'
time = datetime.now().strftime('%Y-%m-%d %H:%M')
file = open(f'/scratch/6256481/results/all/{c}_{i}_{time}.txt', 'w')
file.write( f"FAILED\n{code_failure}")
file.close()