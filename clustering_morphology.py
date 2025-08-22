
# import packages
import xarray as xr
import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools

from datetime import datetime as dt
from glob import glob
from tqdm import tqdm

# import partial packages
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, davies_bouldin_score
from sklearn.utils import resample
from sklearn.model_selection import train_test_split, ParameterGrid
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture

from functools import reduce

from hdbscan import HDBSCAN
from hdbscan.validity import validity_index
from multiprocessing import Pool
from functools import partial

directory = '/scratch/6256481/'

# import sys
# sys.path.insert(0, directory + f'python/py_code/')

###############################################################
# Functions
###############################################################
def scale(df, scaleCols):
    dfCluster = df[scaleCols].copy()
    dfCluster = dfCluster.dropna()
    X = dfCluster[scaleCols]
    # === Standardize features ===
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    return X_scaled

def resampleSilhouette(algo, X, sampleSize,random_state =42, printScore = False):
    labels             = algo.fit_predict(X)
    X_sample, y_sample = resample(X, labels, n_samples=sampleSize, random_state=random_state)
    silScore           = silhouette_score(X_sample, y_sample)
    davScore           = davies_bouldin_score(X_sample, y_sample)

    if printScore == True:
        print(f"Silhouette Score (sample size-{sampleSize}, random_state-{random_state}):", np.round(silScore,3), np.round(davScore, 3))
    return silScore, davScore

def stratified_sample(X, labels, n, rs = 42):
    if len(np.unique(labels)) < 2 or all(labels == -1):
        return None, None
    print('stratified_sample', n)
    X_sample, _, y_sample, _ = train_test_split(X, labels, train_size=n, stratify=labels, random_state=rs)
    return X_sample, y_sample

def HDBSCAN_tune(vars, X):
    metric, algo , minP, sampleSize, mcs, ms = vars

    print(f'min_cluster_size: {mcs}, min_samples: {ms}')
    paramsDict = {'min_cluster_size': mcs, 'min_samples': ms}

    if metric == 'minkowski':
        model = HDBSCAN(**paramsDict, metric = metric, gen_min_span_tree=True, algorithm = algo, metric_params={'P': minP})
    else:
        model = HDBSCAN(**paramsDict, metric = metric, gen_min_span_tree=True, algorithm = algo)
    
    fitted_model = model.fit(X)
    labels = fitted_model.labels_
    if len(set(labels)) > sampleSize:
        sampleSize = int(len(set(labels))*1.2)

    if len(set(labels)) > 1 and np.sum(labels != -1) > 10:
        Xs, Ls    = stratified_sample(X, labels, sampleSize)
        Xs        = Xs.values
        try:
            dbcvScore = validity_index(Xs, Ls)
            relvScore = model.relative_validity_
            print(f"\tDBCV: {dbcvScore:.3f}, relative validity: {relvScore}, Num Clusters: {len(set(labels))}, Noise: {np.sum(labels != -1)}")

        except:
            dbcvScore = relvScore = np.nan
    else:
        dbcvScore = relvScore = np.nan


    return metric, mcs, ms, dbcvScore, relvScore, len(set(labels)), np.sum(labels != -1)

def GMM_tune(vars, X):
        init_type, max_iter, covType, init, cluster=  vars

        gmm = GaussianMixture(n_components  = cluster, 
                            covariance_type = covType,
                            max_iter        = max_iter,
                            init_params     = init_type,
                            random_state    = 42,
                            n_init          = init)
        try:
            gmm.fit(X)

            labels = gmm.predict(X)
            # 4. Calculate BIC score (lower is better)
            bic = gmm.bic(X)
            aic_score = gmm.aic(X)
            
            s,d = [], []
            for rs in random_states:


                _, labelCount = np.unique(labels, return_counts=True)
                if all(labelCount>1):
                    Xs, Ls = stratified_sample(X, labels, sampleSize, rs)
                    # 5. Calculate Silhouette score (range -1 to 1, higher is better)
                    s.append(silhouette_score(Xs, Ls))
                    d.append(davies_bouldin_score(Xs, Ls))
                else:
                    s.append(np.nan), d.append(np.nan)
            sil_score = np.nanmean(s)
            dav_score = np.nanmean(d)
            # print(f"BIC score: {bic:.2f}")
            # print(f"AIC: {aic_score:.3f}")
            print(f"Silhouette Score: {sil_score:.3f}")
            print(f"Davies bouldin Score: {dav_score:.3f}")
        except:
            bic = aic_score = sil_score = dav_score = np.nan
        return init_type, max_iter, covType, init, cluster, bic, aic_score, sil_score, dav_score

###############################################################
# Open files
###############################################################
f = directory + f'results/single_values/global_50.nc'
ds = xr.open_dataset(f)

subsetCols = ['sin', 'bendSin', 'ang', 'apex', 'bendAng','bendWidths', 'bendWidthRatio', 'widthRatio', 'combined_reach_id', 'file_cont', 'file_id']
dsSubset = ds[subsetCols]
df = dsSubset.to_dataframe().reset_index()

df['bendAng'] = abs(df['bendAng'])

df['dimApex'] = df['apex'] / df['bendWidths']
df['bendWidthRatio'] = 1 / df['bendWidthRatio']
df = df[(df['bendSin'] < 6) & (df['dimApex'] < 35) & (df['bendAng'] < 0.1) & (df['bendWidthRatio'] > 0.02)]

# df = df.iloc[0:4000]


dfG = df.groupby(['file_cont', 'combined_reach_id'], as_index = False).agg({
    'sin'       : 'first' ,
    'widthRatio': 'first' ,
    'dimApex'   : 'mean'  ,
    'ang'       : 'mean'  })
FMPCols         = [['sin', 'ang', 'widthRatio', 'dimApex'], ['bendSin', 'bendAng', 'bendWidthRatio', 'dimApex']]
clusterDataType = ['reach', 'bend']

for i, dfT in enumerate([dfG, df]):
    clusters      = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    sampleSize    = 400
    random_states = [20,43,50]
    print()
    print()
    print()
    print()
    print(f'Confinement Height: {FMPCols[i]}')
    print()

    ###############################################################
    # Select columns
    ###############################################################
    for p in [True, False]:
        clusterCols = FMPCols[i] # select grouped clusterCols or bend values
        dfCluster   = dfT.copy().dropna(subset = clusterCols)
        ###############################################################
        # Scale selected columns
        ###############################################################
        XScaled  = scale(dfT, clusterCols)
        
        ###############################################################
        # PCA
        ###############################################################
        print(f'PCA {p}')
        if p == True:
            numPCA = 4
            pca = PCA(n_components=numPCA)  # e.g., reduce to 10 components
            X_pca = pca.fit_transform(XScaled)
            X_pca = X_pca[:,:numPCA]

            clusterCols            = [f'PCA{P}' for P in range(1,numPCA+1)]
            dfCluster[clusterCols] = XScaled[:, :numPCA]

        ###############################################################
        # Kmeans
        # ###############################################################
        print('start Kmeans')
        silScore, davScore = [], []
        
        scores = []
        for j,c in enumerate(clusters):
            print(c)
            sL, dL = [], []
            for rs in random_states:
                kmeans = KMeans(n_clusters=c,init='k-means++',random_state=rs)

                s, d = resampleSilhouette(kmeans, dfCluster[clusterCols]  , sampleSize, rs,False)
                sL.append(s)
                dL.append(d)
            silScore.append(np.mean(sL))
            davScore.append(np.mean(dL))
            scores.append([c, np.mean(sL), np.mean(dL)])
            print(f"Silhouette Score (sample size-{sampleSize}, random_state-{rs}):", np.round(silScore[j],3), np.round(davScore[j], 3))
            print()
        dfKMEANSScores = pd.DataFrame(scores, columns = ['clusters', 'Silhouette', 'DaviesB'])
        dfKMEANSScores.to_csv(directory + f'results/KmeanScores_morphology_{clusterDataType[i]}_{p}.csv')


        ###############################################################
        # PCA GMM
        ###############################################################
        print('GMM')
        init_type = ['kmeans']
        max_iter  = [100]
        covType   = ['diag', 'tied']
        init      = [5]
        clusters  = [1,2,3]
        runs = list(itertools.product(*[max_iter, covType, init, init_type, clusters]))

        if __name__ == '__main__':
            print('__main__')
            func = partial(GMM_tune, X=dfCluster[clusterCols])  # bind df once
            with Pool(10) as pl:
                res = pl.imap(func, runs)
                pl.close()
                pl.join()

        dfGMMScores = pd.DataFrame(list(res), columns = [ 'init_type', 'max_iter', 'covType', 'init', 'cluster','bic', 'aic_score', 'sil_score', 'dav_score'])
        dfGMMScores.to_csv(directory + f'results/GMMSCORES_morphology_{clusterDataType[i]}_{p}.csv')
        print('Done GMM')
        ###############################################################
        # PCA HDBSCAN
        ###############################################################
        print('HBDSCAN')
        dt1 = dt.now()
        runs1 = list(itertools.product(*[['chebyshev'], ['best'], [0], [sampleSize], 
                                        [100, 120, 140, 160,180,200,220,240],
                                        [80, 100, 120, 140, 160,180,200,220]]))
        runs2 = list(itertools.product(*[['manhattan'], ['best'], [0], [sampleSize], 
                                        [60,70,80],
                                        [50,60,70,80]]))
        runs = runs1+runs2
        runs = [r for r in runs if r[-1] <= r[-2]]
        print(runs)
        if __name__ == '__main__':
            print('__main__')
            func = partial(HDBSCAN_tune, X=dfCluster[clusterCols])  # bind df once
            with Pool(10) as pl:
                res = pl.imap(func, runs)
                pl.close()
                pl.join()


    dfHBSCANScores = pd.DataFrame(list(res), columns = ['metric','min_cluster_size', 'min_samples', 'DBCV', 'relVal', 'num_clusters', 'noise'])
    dfHBSCANScores.to_csv(directory + f'results/HDBSCAN_morphology_{clusterDataType[i]}_{p}.csv')
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