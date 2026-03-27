
# import packages
import xarray as xr
import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
from sklearn.mixture import GaussianMixture, BayesianGaussianMixture
# from sklearn.mixture import 

from functools import reduce, partial
import itertools
from multiprocessing import Pool

from joblib import Parallel, delayed

from calc_functions import scale

directory = '/scratch/6256481/'

from clustering import KmeansTune, GMM_tune#, HDBSCAN_tune

###############################################################
# Formulas
###############################################################

def open_dataset_confinement_clustering(ch):
    f = directory + f'results/single_smoothed/global_50_0{ch}_smoothed.nc'
    ds = xr.open_dataset(f)
    df = ds.to_dataframe().reset_index()
    return df

###############################################################
# Tuning for different Conf heights
###############################################################
if __name__ == '__main__':
    clusters      = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    sampleSize    = 40000
    random_states = [20,43,50]
    for ch in [2,3,4]:
        print(ch)

        print()
        print()
        print()
        print()
        print(f'Confinement Height: {ch}')
        print()
        scoresKmeans, scoresGMM, scoresHDB = [], [], []

        df = open_dataset_confinement_clustering(ch)
        ###############################################################
        # Select columns
        ###############################################################
        slopeCols    = ['slope_out_normalized', 'slope_right_normalized', 'slope_left_normalized', 'slope_inn_normalized']
        slopeSCols   = ['slope_out_normalized_smooth', 'slope_right_normalized_smooth', 'slope_left_normalized_smooth', 'slope_inn_normalized_smooth']
        slopeSTCols  = ['slope_out_normalized_smoothSTD', 'slope_right_normalized_smoothSTD', 'slope_left_normalized_smoothSTD', 'slope_inn_normalized_smoothSTD']

        ###############################################################
        # Scale selected columns
        ###############################################################
        # X_scaledO   = scale(df, slopeCols)
        # X_scaledS   = scale(df, slopeSCols)
        # X_scaledST  = scale(df, slopeSCols + slopeSTCols)
        X_scaledSO  = scale(df, slopeCols + slopeSCols)
        # X_scaledSOT = scale(df, slopeCols + slopeSCols + slopeSTCols)
        ###############################################################
        # PCA
        ###############################################################
        pca = PCA(n_components=8)  # e.g., reduce to 10 components
        X_pca = pca.fit_transform(X_scaledSO)
        X_pca = X_pca[:,:3]

        # dfCluster = df[slopeCols + slopeSCols + slopeSTCols]
        dfCluster = df.copy().dropna(subset = slopeCols + slopeSCols)
        dfCluster[['PCA1', 'PCA2', 'PCA3']] = X_pca[:, :3]
        clusterCols = dfCluster.columns[dfCluster.columns.str.startswith('PCA')]

        # for alg in ['tied', 'full', 'diag']:
        #     # test for number of clusters
        #     model = BayesianGaussianMixture(
        #         n_components=40,  # Upper bound
        #         covariance_type=alg,
        #         weight_concentration_prior_type='dirichlet_process',
        #         weight_concentration_prior=1e-2,
        #         max_iter=1000,
        #         random_state=42
        #         )

        #     print(model.fit(dfCluster[clusterCols]))
        #     # Weights represent the proportion of samples in each component
        #     active_components = np.sum(model.weights_ > 1e-2)
        #     print(f"Estimated number of clusters: {active_components}")
        #     print()
        # print()
        # print()

        ###############################################################
        # Kmeans
        # ###############################################################
        print('start Kmeans')           
        S = KmeansTune(dfCluster, clusterCols, clusters, sampleSize, random_states)
        S = [s +comb for s in S]
        scoresKmeans.extend(S)
        # dfKMEANSScores = pd.DataFrame(scores, columns = ['clusters', 'Silhouette', 'DaviesB'])
        # dfKMEANSScores.to_csv(directory + f'results/KmeanScores_confinement_0{ch}.csv')




        ###############################################################
        # PCA GMM
        ###############################################################
        # print('GMM')
        # def GMM_tune(X, n_components, n_init, covType, mi, ip):

        #     gmm = GaussianMixture(n_components=n_components, 
        #                         covariance_type=covType,
        #                         max_iter = mi,
        #                         init_params=ip,
        #                         random_state=42, n_init = n_init)
        #     gmm.fit(X)

        #     labels = gmm.predict(X)
        #     # 4. Calculate BIC score (lower is better)
        #     bic = gmm.bic(X)
        #     aic_score = gmm.aic(X)
            
        #     s,d = [], []
        #     for rs in random_states:
        #         Xs, Ls = stratified_sample(X, labels, sampleSize, rs)
        #         # 5. Calculate Silhouette score (range -1 to 1, higher is better)
        #         s.append(silhouette_score(Xs, Ls))
        #         d.append(davies_bouldin_score(Xs, Ls))
        #     sil_score = np.mean(s)
        #     dav_score = np.mean(d)
        #     # print(f"BIC score: {bic:.2f}")
        #     # print(f"AIC: {aic_score:.3f}")
        #     print(f"Silhouette Score: {sil_score:.3f}")
        #     print(f"Davies bouldin Score: {dav_score:.3f}")
        #     return bic, aic_score, sil_score, dav_score



        scores = []
        for init_p in ['k-means++', 'kmeans']:
            for max_iter in [200, 500]:
                for covType in ['tied', 'diag']:
                    for init in [5, 15, 25, 50]:
                        for n in clusters:
                        
                            b, a, s, d = GMM_tune(dfCluster[clusterCols], n, init, covType, max_iter, init_p)
                            scores.append([max_iter,covType, init, init_p, n,b,a,s,d])
                            # print(max_iter,covType, init, n)
                            # time.sleep(2)

        dfGMMScores = pd.DataFrame(scores, columns = ['max_iter', 'cov_type', 'init', 'init_p', 'clusters', 'BIC', 'AIC', 'Silhouette', 'DaviesB'])
        dfGMMScores.to_csv(directory + f'results/GMMSCORES_confinement_0{ch}.csv')
        print('Done GMM')


        print('GMM')
        init_type = ['k-means++', 'kmeans']
        max_iter  = [200,500]
        covType   = ['diag', 'tied']
        init      = [5, 15, 25, 50]
        runs = list(itertools.product(*[max_iter, covType, init, init_type, clusters]))

        if __name__ == '__main__':
            print('__main__')
            func = partial(GMM_tune, X=dfCluster[clusterCols])  # bind df once
            with Pool(10) as pl:
                res = pl.imap(func, runs)
                pl.close()
                pl.join()

        dfGMMScores = pd.DataFrame(list(res), columns = [ 'init_type', 'max_iter', 'covType', 'init', 'cluster','bic', 'aic_score', 'sil_score', 'dav_score'])
        dfGMMScores.to_csv(directory + f'results/GMMSCORES_confinement_0{ch}.csv')
        print('Done GMM')


        ###############################################################
        # PCA HDBSCAN
        ###############################################################
        # print('HBDSCAN')
        # dt1 = dt.now()

        # print('chebyshev')
        # HDBSCAN_tune_PCA(dfCluster[clusterCols], 'chebyshev', 'best')

        # print('canberra')
        # HDBSCAN_tune_PCA(dfCluster[clusterCols], 'canberra', 'best')

        # print('minkowski 1')
        # HDBSCAN_tune_PCA(dfCluster[clusterCols], 'minkowski', 'best', 1)
        # print('minkowski 2')
        # HDBSCAN_tune_PCA(dfCluster[clusterCols], 'minkowski', 'best', 2)
        # print('minkowski 3')
        # HDBSCAN_tune_PCA(dfCluster[clusterCols], 'minkowski', 'best', 3)

        # print('euclidean')
        # HDBSCAN_tune_PCA(dfCluster[clusterCols], 'euclidean', 'best')
        # dt4 = dt.now()

        # print('manhattan')
        # HDBSCAN_tune_PCA(dfCluster[clusterCols], 'manhattan', 'best')


        # all need algorithm generic and this causes memory issues
        # Cosine creates pairwise datamatrix and this becomes to large
        # correlation creates pairwise datamatrix and this becomes to large
        # seuclidean
        # mahalanobis
        # callable

        # haversine --> needs geographic


    ##########################################################################################
    ##########################################################################################
    # Merged clustering CH 2-3-4
    ##########################################################################################
    ##########################################################################################
    # print('Merged Confinement height approach')
    # def scale(df, scaleCols):
    #     X = df[scaleCols]
    #     # === Standardize features ===
    #     scaler = StandardScaler()
    #     X_scaled = scaler.fit_transform(X)
    #     return X_scaled

    # def open_dataset(CH):
    #     # Open netcdf file
    #     f = directory + f'results/single_smoothed/global_50_0{CH}_smoothed.nc'
    #     ds = xr.open_dataset(f)

    #     # transform to dataframe
    #     df = ds.to_dataframe().reset_index()

    #     # define clustering columns
    #     slopeCols   = ['slope_out_normalized', 'slope_right_normalized', 'slope_left_normalized', 'slope_inn_normalized']
    #     slopeSCols  = ['slope_out_normalized_smooth', 'slope_right_normalized_smooth', 'slope_left_normalized_smooth', 'slope_inn_normalized_smooth']

    #     # drop rows with na values in any of the clustercolumns
    #     dfCluster = df.copy().dropna(subset = slopeCols + slopeSCols )
    #     dfCluster = dfCluster.rename(columns = {col: col + f'_0{CH}' for col in slopeCols + slopeSCols})

    #     return dfCluster

    # def apply_GMM(df, cols, numClusters, CH):

    #     # GMM settings
    #     max_iter = 100
    #     n_init = 5
    #     covariance_type = 'tied'
    #     # select model and fit data
    #     gmm = GaussianMixture(n_components=numClusters, covariance_type=covariance_type, random_state=42, n_init = n_init, max_iter = max_iter)
    #     gmm.fit(df[cols])
    #     print('done Fitting')
    #     labelsHard = gmm.predict(df[cols])
    #     labelsSoft = gmm.predict_proba(df[cols])

    #     df[f'GMM_pca_{numClusters}_hard_0{CH}']=  labelsHard

    #     for i in range(numClusters):
    #         df[f'GMM_pca_{numClusters}_{i}_soft_0{CH}']=  labelsSoft[:,i]
    #     return df

    # dfL = []
    # for ch in [2,3,4]:
    #     df = open_dataset(ch)
    #     dfL.append(df)

    # for i, ch in enumerate([2,3,4]):
    #         cols = ['bendID', 'file', 
    #                 f'slope_out_normalized_0{ch}'  , f'slope_inn_normalized_0{ch}'  , 
    #                 f'slope_left_normalized_0{ch}' , f'slope_right_normalized_0{ch}',
    #                 f'slope_out_normalized_smooth_0{ch}'  , f'slope_inn_normalized_smooth_0{ch}'  , 
    #                 f'slope_left_normalized_smooth_0{ch}' , f'slope_right_normalized_smooth_0{ch}']

    #         dfL[i] = dfL[i][cols]
    # df = reduce(lambda left, right: pd.merge(left, right, on=['file', 'bendID'], how='inner'), dfL)    


    # clusters      = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20, 21, 22, 23,24,25,26,27,28,29,30]
    # sampleSize    = 40000
    # random_states = [20,43,50]


    # ###############################################################
    # # Select columns
    # ###############################################################
    # F = ['slope_out_normalized', 'slope_inn_normalized','slope_left_normalized', 'slope_right_normalized']
    # cols = []
    # for ch in [2,3,4]:
    #     cols.extend([f + f'_0{ch}' for f in F])
    #     cols.extend([f + f'_smooth_0{ch}' for f in F])

    # ###############################################################
    # # Scale selected columns
    # ###############################################################
    # X_scaled  = scale(df, cols)

    # ###############################################################
    # # PCA
    # ###############################################################
    # print('PCA')
    # numPCA = 5
    # pca = PCA(n_components=20)  # e.g., reduce to 10 components
    # X_pca = pca.fit_transform(X_scaled)
    # X_pca = X_pca[:,:numPCA]


    # PCACols = [f'PCA{P}' for P in range(1,numPCA+1)]
    # dfCluster = df.copy().dropna(subset = cols)
    # dfCluster[PCACols] = X_pca[:, :numPCA]

    # print(PCACols)
    # ###############################################################
    # # Kmeans
    # # ###############################################################
    # print('start Kmeans')
    # silScore, davScore, scores = [], [], []
    # for i,c in enumerate(clusters):
    #     print(c)
    #     sL, dL = [], []
    #     for rs in random_states:
    #         kmeans = KMeans(n_clusters=c,init='k-means++',random_state=rs)

    #         s, d = resampleSilhouette(kmeans, dfCluster[PCACols]  , sampleSize, rs,False)
    #         sL.append(s)
    #         dL.append(d)
    #     silScore.append(np.mean(sL))
    #     davScore.append(np.mean(dL))
    #     scores.append([c, np.mean(sL), np.mean(dL)])
    #     print(f"Silhouette Score (sample size-{sampleSize}, random_state-{rs}):", np.round(silScore[i],3), np.round(davScore[i], 3))
    #     print()
    # dfKMEANSScores = pd.DataFrame(scores, columns = ['clusters', 'Silhouette', 'DaviesB'])
    # dfKMEANSScores.to_csv(directory + f'results/KmeanScores_confinement_merged.csv')


    # ###############################################################
    # # PCA GMM
    # ###############################################################
    # print('GMM')
    # def GMM_tune(X, n_components, n_init, covType, mi):
    #     print(n_components, n_init, covType, mi)
    #     gmm = GaussianMixture(n_components=n_components, 
    #                         covariance_type=covType,
    #                         max_iter = mi,
    #                         init_params='kmeans',
    #                         random_state=42, n_init = n_init)
    #     gmm.fit(X)

    #     labels = gmm.predict(X)
    #     # 4. Calculate BIC score (lower is better)
    #     bic = gmm.bic(X)
    #     aic_score = gmm.aic(X)
        
    #     s,d = [], []
    #     for rs in random_states:
    #         Xs, Ls = stratified_sample(X, labels, sampleSize, rs)
    #         # 5. Calculate Silhouette score (range -1 to 1, higher is better)
    #         s.append(silhouette_score(Xs, Ls))
    #         d.append(davies_bouldin_score(Xs, Ls))
    #     sil_score = np.mean(s)
    #     dav_score = np.mean(d)
    #     # print(f"BIC score: {bic:.2f}")
    #     # print(f"AIC: {aic_score:.3f}")
    #     print(f"Silhouette Score: {sil_score:.3f}")
    #     print(f"Davies bouldin Score: {dav_score:.3f}")
    #     return bic, aic_score, sil_score, dav_score
    #     # dbcv_score = validity_index(Xs, Ls, metric='euclidean')
    #     # print(f"DBCV score: {dbcv_score:.3f}")


    # scores = []
    # for max_iter in [100]:
    #     for covType in ['diag', 'tied']:
    #         for init in [5]:
    #             for n in clusters:
                
    #                 b, a, s, d = GMM_tune(dfCluster[PCACols], n, init, covType, max_iter)
    #                 scores.append([max_iter,covType, init, n,b,a,s,d])
    #                 # print(max_iter,covType, init, n)
    #                 # time.sleep(2)

    # dfGMMScores = pd.DataFrame(scores, columns = ['max_iter', 'cov_type', 'init', 'clusters', 'BIC', 'AIC', 'Silhouette', 'DaviesB'])
    # dfGMMScores.to_csv(directory + f'results/GMMSCORES_confinement_merged.csv')
    # print('Done GMM')



