from run_confinement_values import create_apex_val_dataframe
from glob import glob
from multiprocessing import Pool
import numpy as np
import os

directory = '/scratch/6256481/'
confFactor = 50
files = np.sort(glob(directory + f'results/all/??_??_{confFactor}.csv'))

removeFiles = glob(directory + f'results/single_values/??_??_{confFactor}.csv')
for rmf in removeFiles:
    os.remove(rmf)

# single test
# print(files[4])
# create_apex_val_dataframe(files[4])


if __name__ == '__main__':
    with Pool(10) as p:
        p.imap(create_apex_val_dataframe, files)
        p.close()
        p.join()

print('Open_to_single_apex Finished')