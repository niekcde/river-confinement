if __package__ in (None, ""):
    import os
    import sys

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "pipeline"

################################################
# Import modules
################################################
import numpy as np
from glob import glob
from datetime import datetime as dt

from .build_step3_single_values import run_step3_single_values
from .build_step4_confinement_factor import build_step4_confinement_factor
from .build_step5_confinement_outputs import run_step5_confinement_outputs
from .build_step6_aggregates import run_step6_aggregates

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
    step3_outputs = run_step3_single_values(
        workers=10,
        conf_factor=crossFactor,
    )
    
    dt2 = dt.now()
    print(f'Open_to_single_apex Finished: {dt2-dt1}')
    allResultFiles = np.array([str(path) for path in step3_outputs])

################################################
# get confinement factor values
################################################
if createNewFactor == True:
    print('Start calc confinement factor')
    build_step4_confinement_factor(
        input_files=allResultFiles,
        output_path=directory + f'results/reference_tables/confinement_factor_{crossFactor}.csv',
        y1=conFactor[0],
        y2=conFactor[1],
    )
    
    dt3 = dt.now()
    dtm = dt2 if createNewSingleVal == True else dt1
    print(f'Open_to_single_apex Finished: {dt3 - dtm}')

################################################
# get confinement values
################################################
print('Start calc confinement values')
for hf in heightFactor:
    run_step5_confinement_outputs(
        input_files=allResultFiles,
        workers=10,
        conf_factor=crossFactor,
        height_factor=hf,
        directory=directory,
    )

    run_step6_aggregates(
        cross_factor=crossFactor,
        height_factor=hf,
        directory=directory,
    )

    print(f'run confinement heightfactor {hf} finished')
print(f'run confinement finished Finished: {dt.now() - dt1}')
