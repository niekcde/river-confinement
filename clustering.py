import numpy as np
import pandas as pd

from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score, davies_bouldin_score
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA


# from hdbscan import HDBSCAN
# from hdbscan.validity import validity_index

from calc_functions import log_var, scale, stratified_sample


def HDBSCAN_tune(vars, X):
    metric, algo , minP, sampleSize, mcs, ms = vars
    paramsDict = {'min_cluster_size': mcs, 'min_samples': ms}
    if metric == 'minkowski':
        model = HDBSCAN(**paramsDict, metric = metric, gen_min_span_tree=True, algorithm = algo, metric_params={'P': minP})
    else:
        model = HDBSCAN(**paramsDict, metric = metric, gen_min_span_tree=True, algorithm = algo)
    
    fitted_model = model.fit(X)
    labels       = fitted_model.labels_
    if len(set(labels)) > sampleSize:
        sampleSize = int(len(set(labels))*1.2)

    if len(set(labels)) > 1 and np.sum(labels != -1) > 10:
        Xs, Ls    = stratified_sample(X, labels, 1, sampleSize)
        Xs        = Xs.values
        try:
            dbcvScore = validity_index(Xs, Ls)
            relvScore = model.relative_validity_
            # print(f"\tDBCV: {dbcvScore:.3f}, relative validity: {relvScore}, Num Clusters: {len(set(labels))}, Noise: {np.sum(labels != -1)}")

        except:
            dbcvScore = relvScore = np.nan
    else:

        dbcvScore = relvScore = np.nan


    return metric, mcs, ms, dbcvScore, relvScore, len(set(labels)), np.sum(labels != -1)

def GMM_tune(vars, X):
    init_type, max_iter, covType, init, cluster, sampleSize, random_states =  vars
    
    gmm = GaussianMixture(n_components  = cluster, 
                            covariance_type = covType,
                            max_iter        = max_iter,
                            init_params     = init_type,
                            random_state    = 42,
                            n_init          = init)
    gmm.fit(X)
    labels = gmm.predict(X)
    bic = gmm.bic(X)
    aic_score = gmm.aic(X)
    
    s,d = [], []
    for rs in random_states:
        _, labelCount = np.unique(labels, return_counts=True)
        if all(labelCount>1):
            Xs, Ls = stratified_sample(X, labels,cluster, sampleSize, rs)
            if Xs is None:
                s.append(np.nan)
                d.append(np.nan)
            else:
                s.append(silhouette_score(Xs, Ls))
                d.append(davies_bouldin_score(Xs, Ls))
        else:
            s.append(np.nan), d.append(np.nan)
    sil_score = np.nanmean(s)
    dav_score = np.nanmean(d)

    return init_type, max_iter, covType, init, cluster, bic, aic_score, sil_score, dav_score

def applyGMM(df, cols, cluster, init_type, max_iter, covType, init):

    gmm = GaussianMixture(n_components    = cluster, 
                          covariance_type = covType, 
                          random_state    = 42, 
                          n_init          = init,
                          init_params     = init_type, 
                          max_iter        = max_iter)
    gmm.fit(df[cols])

    labelsHard = gmm.predict(df[cols])
    labelsSoft = gmm.predict_proba(df[cols])

    df[f'GMM_{cluster}_hard'] = labelsHard

    for c in range(cluster):
        df[f'GMM_{cluster}_{c}_soft']=  labelsSoft[:,c]
    return df

def applyKmeans(df, cols, cluster, init = 'k-means++'):
    df = df.copy()

    kmeans = KMeans(n_clusters=cluster,init=init,random_state=42)
    labels = kmeans.fit_predict(df[cols])
    labels = np.array(labels)
    df[f'kmeans_{cluster}'] = labels
    return df


def KmeansTune(df, cols, clusters, sampleSize, randomStates):
    scores = []
    for ci, c in enumerate(clusters):
        silScore, davScore = [],[]
        for rs in randomStates:
            kmeans = KMeans(n_clusters=c,init='k-means++',random_state=rs)
            labels = kmeans.fit_predict(df[cols])
            Xs, Ys = stratified_sample(df[cols], labels, c, sampleSize, rs = rs)
            silScore.append(silhouette_score(Xs, Ys))
            davScore.append(davies_bouldin_score(Xs, Ys))
        scores.append([c, np.mean(silScore), np.mean(davScore)])
    return scores


# Suppose you have original feature names
def print_pca_loadings(pca, colsNames):
    # Create a DataFrame of loadings
    loadings = pd.DataFrame(
        data=pca.components_.T,
        index=colsNames,
        columns=[f"PC{i+1}" for i in range(pca.n_components_)]
    )
    print(loadings.round(3))


def data_prep_clustering(df,cols,applyCPA, numPCA, unknownPCA):
    dfCluster = df.copy().dropna(subset = cols)
    X_scaled = scale(dfCluster, cols)
    if applyCPA == True:
        pca = PCA(n_components=len(cols))  # e.g., reduce to 10 components
        Xpca = pca.fit_transform(X_scaled)
        
        if unknownPCA == True:
            print(pca.explained_variance_ratio_)
            numPCA = len(cols)
            pcaCols = [f'PCA{P}' for P in range(1,numPCA+1)]
            dfCluster[pcaCols] = Xpca[:, :numPCA]
        else:
            pcaCols = [f'PCA{P}' for P in range(1,numPCA+1)]
            dfCluster[pcaCols] = Xpca[:, :numPCA]
        return dfCluster, pca, Xpca
    else:
        cols = [f'{c}_scaled' for c in cols]
        dfCluster[cols] = X_scaled
        return dfCluster