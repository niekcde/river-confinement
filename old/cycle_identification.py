import numpy as np
import geopandas as gpd
import pandas as pd
import shapely

import networkx as nx

from itertools import combinations
# import tqdm



def get_nodes(DF):
    df = DF.copy()
    Nodes = gpd.GeoDataFrame(columns = ['id', 'reach_id', 'geometry'])

    # for index1, row1 in tqdm.tqdm(df.iterrows(), total=df.shape[0]):
    for index1, row1 in df.iterrows():

        
        dists = row1.geometry.distance(df.geometry).values

        df.loc[:, 'dist'] = dists
        df = df.sort_values('dist')
        

        potential_edges = df[df.dist < 2e2]
        if potential_edges.shape[0] == 1:
            continue;

 
        pot = potential_edges.copy()

        E = pot[pot.reach_id != row1.reach_id]
        I = 0
        

        for i,e in E.iterrows():


            Ndist = e.geometry.distance(E.geometry)
            Sel = E[Ndist < 2e2]

            P = shapely.ops.nearest_points(e.geometry, row1.geometry)[0]

            if np.isnan(Nodes.id.max()):
                Nid = [1]
            else:
                Nid = [Nodes.id.max() + 1]

            Rid = np.append(Sel.reach_id.values, row1.reach_id)

            Nid = Nid*len(Rid)
            P   = [P]*len(Rid)

            if Sel.shape[0] == 1:

                N = gpd.GeoDataFrame({'id': Nid,'reach_id':Rid, 'geometry': P})
            else:
                N = gpd.GeoDataFrame({'id': Nid,'reach_id':Rid, 'geometry': P})

            Nodes = pd.concat([Nodes, N], ignore_index = True)

            E = E.drop(index = Sel.index)
            if E.shape[0] == 0:
                break;
            

    return Nodes


def combine_remove_Nodes(nodes):
    Nodes = nodes.copy()
    Nids  = Nodes.id.unique()
    Nodes.loc['Nid', :] = pd.Series(dtype='int')

    for N in range(len(Nids)):
        Nid = Nids[N]
        
        TN = Nodes[Nodes.id == Nid].iloc[0]

        potNE = TN.geometry.distance(Nodes.geometry)
        potNE = Nodes[potNE < 2e2]
        Nodes.loc[potNE.index, 'Nid'] = potNE.id.min()

    Nodes = Nodes.drop_duplicates(['reach_id', 'Nid'])
    Nodes = Nodes.dropna()
    return Nodes


def get_edges(nodes, df):
    Nodes = nodes.copy()
    edges = []
    ids   = []
    for r in Nodes.reach_id.unique():
        rs    = Nodes[Nodes.reach_id == r]
        reach = df[df.reach_id == r]
        # Generate combinations to check if reach has node half way
        combinations_list = list(combinations(rs.Nid, 2))

        if len(combinations_list) == 3:
            
            D = 0
            for c in combinations_list:
                D1 = reach.iloc[0].geometry.project(rs[rs.Nid == c[0]].iloc[0].geometry)
                D2 = reach.iloc[0].geometry.project(rs[rs.Nid == c[1]].iloc[0].geometry)
                if (abs(D1 - D2) > D):
                    D  = abs(D1 - D2)
                    dc = c
            combinations_list.remove(dc)
        elif len(combinations_list) > 3:
            
            dists = [] 
            for c in rs.Nid.values:
                dists.append(reach.iloc[0].geometry.project(rs[rs.Nid == c].iloc[0].geometry))
            dists.sort() 
            corr_dists = np.diff(dists)
            rem = []           
            for c in combinations_list:
                D1 = reach.iloc[0].geometry.project(rs[rs.Nid == c[0]].iloc[0].geometry)
                D2 = reach.iloc[0].geometry.project(rs[rs.Nid == c[1]].iloc[0].geometry)
                if abs(D1-D2) not in corr_dists:
                    rem.append(c)
            combinations_list = [X for X in combinations_list if X not in rem]

                             
    

        edges.extend(combinations_list)
        for i in range(len(combinations_list)):
            ids.append(rs.iloc[0].reach_id)
    

    dfE = pd.DataFrame({'reach_id':ids, 'Edge': edges})
    return dfE

def get_cycles(df, plot = False):
    Nodes = get_nodes(df)
    Nodes = combine_remove_Nodes(Nodes)
    edges = get_edges(Nodes, df)
    
    G = nx.Graph()
    G.add_edges_from(edges.Edge)

    cycles = nx.cycle_basis(G)
    rids = []
    cids = []
    for i, c in enumerate(cycles):
        rid = Nodes[Nodes.Nid.isin(c)].reach_id.unique()

        rids.extend(rid)

        cid = [i] *len(rid)
        cids.extend(cid)

    df_cycle = pd.DataFrame({'reach_id':rids, 'Cycle':cids})
    
    col = 'Cycle'
    new_df = pd.DataFrame(df_cycle.groupby('reach_id')[col].agg(lambda x: x.tolist()).tolist()).replace({None: np.nan})
    new_df.columns = [f'{col}{i}' for i in new_df.columns + 1]


    new_df['reach_id'] = df_cycle['reach_id'].drop_duplicates().reset_index(drop=True)



    return new_df