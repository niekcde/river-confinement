if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

################################################
# Import modules
################################################
import pandas as pd
import numpy as np
import os
from glob import glob
from multiprocessing import Pool
from datetime import datetime as dt

from .support import confinement_factor_single_values
from .run_confinement_values import calc_confinement_values, concat_nc_conf_files,\
      concat_reachAveraged, create_apex_val_dataframe

################################################
# Settings
################################################
directory = '/scratch/6256481/'
crossFactor  = 50
heightFactor = [2, 0.5, 1, 1.5, 3, 4, 6, 8, 10, 15]
conFactor    = [50,10]
createNewSingleVal = True
createNewFactor    = True
allResultFiles = np.sort(glob(directory + f'results/single_values/??_??_{crossFactor}.csv'))

dt1 = dt2 = dt3 = dt.now()
print(f'start Code: {dt1}')

################################################
# Create singe bend value files
################################################
if createNewSingleVal == True:
    print('Start transform results in single row values')
    files = np.sort(glob(directory + f'results/all/??_??_{crossFactor}.csv'))
    removeFiles = glob(directory + f'results/single_values/??_??_{crossFactor}.csv')
    for rmf in removeFiles:
        os.remove(rmf)

    if __name__ == '__main__':
        with Pool(10) as p:
            p.imap(create_apex_val_dataframe, files)
            p.close()
            p.join()
    
    dt2 = dt.now()
    print(f'Open_to_single_apex Finished: {dt2-dt1}')
    allResultFiles = np.sort(glob(directory + f'results/single_values/??_??_{crossFactor}.csv'))

################################################
# get confinement factor values
################################################
if createNewFactor == True:
    print('Start calc confinement factor')
    for i, f in enumerate(allResultFiles):
        dfTemp = pd.read_csv(f, dtype = {'include_flag':str, 'calculated':str})
        if i == 0:
            dfA = dfTemp
        else:
            dfA = pd.concat([dfA, dfTemp])

    dfCF = confinement_factor_single_values(dfA, 'bendWidths', conFactor[0], conFactor[1])
    dfCF.to_csv(directory + 'results/confinement_factor.csv')
    
    dt3 = dt.now()
    dtm = dt2 if createNewSingleVal == True else dt1
    print(f'Open_to_single_apex Finished: {dt3 - dtm}')

################################################
# get confinement values
################################################
print('Start calc confinement values')
def run(input):
    file = input[0]
    print(f'run_confinement_values_shell - run: {file[-12:-4]}')
    df = pd.read_csv(file)
    calc_confinement_values(df, file[-12:-4], directory, False, True, crossFactor, input[1])


# dfFiles = pd.read_csv(directory + 'results/reference_tables/file_sorting.csv', index_col = 0)
# sorting = list(dfFiles.sort_values('file')['size'].argsort().values[::-1])
# allResultFiles = np.array(allResultFiles)
# allResultFiles = allResultFiles[sorting]

# print(len(allResultFiles))
# run([allResultFiles[-1], 2])


for hf in heightFactor:
    multiInput = [[arf, hf] for arf in allResultFiles]
    if __name__ == '__main__':
        with Pool(10) as p:
            p.imap(run, multiInput)
            p.close()
            p.join()

    concat_nc_conf_files(directory, crossFactor, hf)
    concat_reachAveraged(directory, crossFactor, hf)

    print(f'run confinement heightfactor {hf} finished')
print(f'run confinement finished Finished: {dt.now() - dt1}')
