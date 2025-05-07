import pandas as pd
import numpy as np

from glob import glob
from multiprocessing import Pool
from datetime import datetime as dt

from support import confinement_factor_single_values
from run_confinement_values import calc_confinement_values, concat_nc_conf_files, concat_reachAveraged



directory = '/scratch/6256481/'
crossFactor = 50
conFactor  = [50,10]
createNewFactor = False
########################
# get confinement factor values
########################
dt1 = dt.now()
print(f'start Code: {dt.now()}')
allResultFiles = np.sort(glob(directory + f'results/single_values/??_??_{crossFactor}.csv'))
if createNewFactor == True:
    for i, f in enumerate(allResultFiles):
        dfTemp = pd.read_csv(f, dtype = {'include_flag':str, 'calculated':str})
        if i == 0:
            dfA = dfTemp
        else:
            dfA = pd.concat([dfA, dfTemp])

    dfCF = confinement_factor_single_values(dfA, 'bendWidths', conFactor[0], conFactor[1])
    dfCF.to_csv(directory + 'results/confinement_factor.csv')

print('start MultiProcess')
def run(file):
    print(f'run_confinement_values_shell - run: {file[-12:-4]}')
    df = pd.read_csv(file)
    calc_confinement_values(df, file[-12:-4], directory, False, True)


print('number of multi process files: ', len(allResultFiles))

# dfFiles = pd.read_csv(directory + 'results/file_sorting.csv', index_col = 0)
# sorting = list(dfFiles.sort_values('file')['size'].argsort().values[::-1])
# allResultFiles = np.array(allResultFiles)
# allResultFiles = allResultFiles[sorting]

# print(allResultFiles[-1])
# run(allResultFiles[-1])

if __name__ == '__main__':
    with Pool(10) as p:
        p.imap(run, allResultFiles)
        p.close()
        p.join()

concat_nc_conf_files(directory, crossFactor)
concat_reachAveraged(directory, crossFactor)

print(f'run confinement finished Finished: {dt.now() - dt1}')