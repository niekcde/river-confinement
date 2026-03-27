
# import packages
import xarray as xr
import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools
from joblib import Parallel, delayed

from datetime import datetime as dt
from glob import glob
from tqdm import tqdm

# import partial packages
from sklearn.decomposition import PCA
from functools import partial

directory = '/scratch/6256481/'
directory = '/Users/6256481/Library/CloudStorage/OneDrive-UniversiteitUtrecht/Desktop/Congres/AGU_2025/presentation/analysis/'
# import functions
from calc_functions import log_var, scale
from clustering import KmeansTune, GMM_tune, HDBSCAN_tune
###############################################################
# Functions
###############################################################
def valid_combination(comb):
    # Extract base variables ignoring "_log"
    base_vars_in_comb = [v.replace('Log','') for v in comb]
    # Check for duplicates → if yes, invalid (same variable both logged/unlogged)
    return len(base_vars_in_comb) == len(set(base_vars_in_comb))

def open_dataset_morphology_clustering():
    # f = directory + f'results/single_values/global_50.nc'
    f = directory + f'global_50.nc'
    ds = xr.open_dataset(f)

    subsetCols = ['sin', 'bendSin', 'ang', 'apex', 'bendAng','bendWidths', 'bendWidthRatio', 'widthRatio', 'combined_reach_id', 'file_cont', 'file_id']
    dsSubset   = ds[subsetCols]
    df         = dsSubset.to_dataframe().reset_index()

    # Transform columns
    df['bendAng']        = abs(df['bendAng'])
    df['dimApex']        = df['apex'] / df['bendWidths']
    df['bendWidthRatio'] = 1 / df['bendWidthRatio']
    # Filter outliers
    df = df[(df['bendSin'] < 6) & 
            (df['dimApex'] < 35) & 
            (df['bendAng'] < 0.1) & 
            (df['bendWidthRatio'] > 0.02)]

    # Create reach averages
    dfG = df.groupby(['file_cont', 'combined_reach_id'], as_index = False).agg({
        'sin'       : 'first' ,
        'widthRatio': 'first' ,
        'dimApex'   : 'mean'  ,
        'bendAng'   : 'mean'  })
    dfG = dfG.rename(columns= {'bendAng':'ang'})

    return df, dfG
###############################################################
# Selection of code only run if actual file is run
###############################################################
if __name__ == "__main__":
    # Open files
    df, dfG         = open_dataset_morphology_clustering()
    FMPCols         = [['sin', 'ang', 'widthRatio', 'dimApex'], ['bendSin', 'bendAng', 'bendWidthRatio', 'dimApex']]
    clusterDataType = ['reach', 'bend']

    clusters      = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    sampleSize    = 40000
    random_states = [20,43,50]

    # dfG           = dfG.iloc[0:1000]
    # clusters      = [3,4]
    # sampleSize    = 990
    # random_states = [20]
    for i, dfT in enumerate([dfG]):
        ###############################################################
        # Select columns
        ###############################################################
        for p in [1]:
            clusterCols = FMPCols[i] # select grouped clusterCols or bend values
            dfCluster   = dfT.copy().dropna(subset = clusterCols)
            
            ###############################################################
            # Log transform selected coloms
            ###############################################################
            dfCluster       = log_var(dfCluster, clusterCols)
            clusterCols_log = [x + 'Log' for x in clusterCols]
            clusterColsAll  = clusterCols + clusterCols_log
            ###############################################################
            # create column combination and start loop
            ###############################################################
            comb_4 = list(itertools.combinations(clusterColsAll, 4))
            valid_comb_4 = [c for c in comb_4 if valid_combination(c)]
            combColNames = [f'col{c}' for c in range(4)]
            scoresKmeans, scoresGMM, scoresHDB = [], [], []
            for comb in tqdm(valid_comb_4):

                comb = list(comb)
                ###############################################################
                # Scale selected columns
                ###############################################################
                XScaled  = scale(dfCluster, comb)
                
                ###############################################################
                # PCA
                ###############################################################
                # print(f'PCA {p}')
                if p in [1,2]:
                    numPCA = 4
                    pca = PCA(n_components=numPCA)  # e.g., reduce to 10 components
                    X_pca = pca.fit_transform(XScaled)
                    X_pca = X_pca[:,:numPCA]

                    clusterCols            = [f'PCA{P}' for P in range(1,numPCA+1)]
                    dfCluster[clusterCols] = X_pca[:, :numPCA]
                else:
                    dfCluster[clusterCols] = XScaled
                ###############################################################
                # Kmeans
                # ###############################################################
                # print('start Kmeans')           
                S = KmeansTune(dfCluster, clusterCols, clusters, sampleSize, random_states)
                S = [s +comb for s in S]
                scoresKmeans.extend(S)

                ###############################################################
                # PCA GMM
                ###############################################################
                init_type = ['kmeans', 'k-means++']
                max_iter  = [100]
                covType   = ['diag', 'tied', 'full']
                init      = [5, 10, 15]
                samplesSL = [sampleSize]
                rsL       = [random_states]

                runs = list(itertools.product(*[init_type, max_iter, covType, init, clusters, samplesSL, rsL]))
                func = partial(GMM_tune, X=dfCluster[clusterCols])  # bind df once

                res = Parallel(n_jobs=5)(delayed(func)(r) for r in runs)
                # print(results)
                GMMResults = [list(x) + comb for x in res]
                scoresGMM.extend(GMMResults)
                ###############################################################
                # HDBSCAN
                ###############################################################
                # print('HBDSCAN')
                runs1 = list(itertools.product(*[['chebyshev'], ['best'], [0], [sampleSize], 
                                                [100, 120, 140, 160,180,200,220,240],
                                                [80, 100, 120, 140, 160,180,200,220]]))
                runs2 = list(itertools.product(*[['manhattan'], ['best'], [0], [sampleSize], 
                                                [60,70,80],
                                                [50,60,70,80]]))
                
                runs = runs1+runs2
                runs = [r for r in runs if r[-1] <= r[-2]]

                func = partial(HDBSCAN_tune, X=dfCluster[clusterCols])  # bind df once

                res = Parallel(n_jobs=5)(delayed(func)(r) for r in runs)

                HDBResults = [list(x) + comb for x in list(res)]
                scoresHDB.extend(HDBResults)

            ###############################################################
            # Save final results
            ###############################################################
            dfKMEANSScores = pd.DataFrame(scoresKmeans, columns = ['clusters', 'Silhouette', 'DaviesB'] + combColNames)
            dfKMEANSScores.to_csv(directory + f'results/tuning/KmeanScores_morphology_{clusterDataType[i]}.csv')

            dfGMMScores = pd.DataFrame(scoresGMM, columns = [ 'init_type', 'max_iter', 'covType', 'init', 'cluster','BIC', 'AIC', 'Silhouette', 'DaviesB']+ combColNames)
            dfGMMScores.to_csv(directory + f'results/tuning/GMMSCORES_morphology_{clusterDataType[i]}.csv')
            
            dfHBSCANScores = pd.DataFrame(scoresHDB, columns = ['metric','min_cluster_size', 'min_samples', 'DBCV', 'relVal', 'num_clusters', 'noise']+ combColNames)
            dfHBSCANScores.to_csv(directory + f'results/tuning/HDBSCAN_morphology_{clusterDataType[i]}.csv')
            
            
            
            
            
            
            
            # print('chebyshev')
            # HDBSCAN_tune(dfCluster[clusterCols], 'chebyshev', 'best', sampleSize=sampleSize)

            # print('canberra')
            # HDBSCAN_tune(dfCluster[clusterCols], 'canberra', 'best', sampleSize = sampleSize)

            # print('minkowski 1')
            # HDBSCAN_tune_PCA(dfCluster[clusterCols], 'minkowski', 'best', 1)
            # print('minkowski 2')
            # HDBSCAN_tune_PCA(dfCluster[clusterCols], 'minkowski', 'best', 2)
            # print('minkowski 3')
            # HDBSCAN_tune_PCA(dfCluster[clusterCols], 'minkowski', 'best', 3)

            # print('euclidean')
            # HDBSCAN_tune(dfCluster[clusterCols], 'euclidean', 'best', sampleSize = sampleSize)
            # dt4 = dt.now()

            # print('manhattan')
            # HDBSCAN_tune(dfCluster[clusterCols], 'manhattan', 'best', sampleSize = sampleSize)


            # all need algorithm generic and this causes memory issues
            # Cosine creates pairwise datamatrix and this becomes to large
            # correlation creates pairwise datamatrix and this becomes to large
            # seuclidean
            # mahalanobis
            # callable

            # haversine --> needs geographic

    print('DONE')