"""Build auditable bend -> combined reach -> SWORD -> HydroBASINS links.

Use the matching submission reach files, not a different SWORD release.
Connected-component networkGraph values are NOT watershed identifiers.
"""
from pathlib import Path
import sqlite3
import json
import pandas as pd
import numpy as np
import geopandas as gpd

TARGETS = {
    'rhine': {'continent':'eu','main_bas':2060023010,'reach_file':'eu_02_reach_new_segments.gpkg'},
    'amazon': {'continent':'sa','main_bas':6060007000,'reach_file':'sa_01_reach_new_segments.gpkg'},
    'ob': {'continent':'as','main_bas':3060001840,'reach_file':'as_00_reach_new_segments.gpkg'},
}

def build_lookup(source, vector_dir, hydro_dir, out):
    source,vector_dir,hydro_dir,out=map(Path,[source,vector_dir,hydro_dir,out])
    out.mkdir(parents=True,exist_ok=True)
    frames=[]
    for cont in ['eu','sa','as']:
        for path in sorted(vector_dir.glob(f'{cont}_*.gpkg')):
            with sqlite3.connect(f'file:{path}?mode=ro',uri=True) as db:
                layer=db.execute('select table_name from gpkg_contents').fetchone()[0]
                d=pd.read_sql_query(f'SELECT combined_reach_id, reach_id, networkGraph, river_name FROM "{layer}"',db)
            d['file']=cont;d['source_file']=path.name;frames.append(d)
    reaches=pd.concat(frames,ignore_index=True)
    reaches['combined_reach_id']=reaches.combined_reach_id.astype('int64')
    reaches['PFAF_ID']=reaches.reach_id.astype('int64')//100000
    hydro=pd.concat([gpd.read_file('zip://'+str((hydro_dir/f'hybas_{c}_lev06_v1c.zip').resolve()),
                                 ignore_geometry=True) for c in ['eu','sa','si']],ignore_index=True)
    hydro=hydro[['PFAF_ID','HYBAS_ID','MAIN_BAS','UP_AREA']]
    assert not hydro.PFAF_ID.duplicated().any()
    joined=reaches.merge(hydro,on='PFAF_ID',how='left',validate='many_to_one')
    joined.MAIN_BAS=joined.MAIN_BAS.astype('Int64')
    grouped=joined.groupby(['file','combined_reach_id'],sort=False)
    combined=grouped.agg(reach_count=('reach_id','size'),mapped_reach_count=('MAIN_BAS','count'),
                         main_bas_count=('MAIN_BAS','nunique'),MAIN_BAS=('MAIN_BAS','first'),
                         source_file=('source_file','first'),networkGraph=('networkGraph','first')).reset_index()
    combined['assignment_status']=np.select(
        [combined.mapped_reach_count.eq(0),combined.main_bas_count.gt(1),combined.mapped_reach_count.lt(combined.reach_count)],
        ['unmapped_basin','multiple_catchments','partly_unmapped'],default='unambiguous')
    combined.loc[combined.assignment_status.ne('unambiguous'),'MAIN_BAS']=pd.NA
    names={v['main_bas']:k for k,v in TARGETS.items()}
    combined['catchment']=combined.MAIN_BAS.map(names)
    assert not combined.duplicated(['file','combined_reach_id']).any()
    combined.to_csv(out/'combined_reach_catchment_lookup.csv.gz',index=False)
    joined.to_csv(out/'sword_reach_hydrobasin_lookup.csv.gz',index=False)
    # Named source reaches verify that each selected MAIN_BAS is the main basin.
    for slug,target in TARGETS.items():
        assert (joined.MAIN_BAS==target['main_bas']).any()
        assert target['reach_file'] in joined.loc[joined.MAIN_BAS==target['main_bas'],'source_file'].unique()
    print('Reading bend IDs from the original GeoPackage',flush=True)
    with sqlite3.connect(f'file:{source}?mode=ro',uri=True) as db:
        bends=pd.read_sql_query('SELECT fid, file, bendID, combined_reach_id_x, combined_reach_id_y FROM global_clusters_bend',db)
    # The file field is the continent; bendID alone is reused across continents.
    bends['combined_reach_id']=bends.combined_reach_id_x.astype('int64')
    prefix=bends.bendID.str.split('_').str[0].astype('int64')
    assert prefix.eq(bends.combined_reach_id).all(), 'bendID prefix disagrees with combined_reach_id_x'
    known=bends.combined_reach_id_y.notna()
    assert bends.loc[known,'combined_reach_id_y'].eq(bends.loc[known,'combined_reach_id_x']).all()
    duplicate_bend_rows = int(bends.duplicated(['file','bendID']).sum())
    linked=bends[['fid','file','bendID','combined_reach_id']].merge(combined,on=['file','combined_reach_id'],how='left',validate='many_to_one')
    linked.assignment_status=linked.assignment_status.fillna('no_matching_combined_reach')
    columns=['fid','file','bendID','combined_reach_id','MAIN_BAS','catchment','assignment_status','source_file']
    # Source rows can repeat a bend key. Their validated combined-reach prefix
    # is identical, so membership is identical; retain one lookup row per key.
    lookup = linked[[c for c in columns if c != 'fid']].drop_duplicates()
    assert not lookup.duplicated(['file','bendID']).any(), 'Conflicting membership for a bend key'
    lookup.to_csv(out/'bend_catchment_lookup.csv.gz',index=False)
    linked.loc[linked.catchment.notna(),columns].to_csv(out/'selected_catchment_bends.csv.gz',index=False)
    ambiguous=joined.merge(combined[['file','combined_reach_id','assignment_status']],on=['file','combined_reach_id'],validate='many_to_one')
    ambiguous=ambiguous[ambiguous.assignment_status.ne('unambiguous') & ambiguous.MAIN_BAS.isin(names)]
    ambiguous.to_csv(out/'ambiguous_target_reaches.csv',index=False)
    info={'source':str(source),'vector_dir':str(vector_dir),'targets':TARGETS,
          'bend_rows':len(bends),'bend_key':['file','bendID'],
          'duplicate_source_bend_key_rows':duplicate_bend_rows,'unique_bend_keys':len(lookup),
          'catchment_bends':linked.catchment.value_counts().to_dict(),
          'catchment_combined_reaches':combined.catchment.value_counts().to_dict(),
          'ambiguous_target_combined_reaches':ambiguous.groupby('file').combined_reach_id.nunique().to_dict(),
          'policy':'Only combined reaches whose constituent SWORD reaches all map to one MAIN_BAS. Ambiguous and incomplete assignments excluded.',
          'hydrobasins_urls':[f'https://data.hydrosheds.org/file/hydrobasins/standard/hybas_{c}_lev06_v1c.zip' for c in ['eu','sa','si']]}
    (out/'catchment_lookup_audit.json').write_text(json.dumps(info,indent=2))
    print(json.dumps(info,indent=2),flush=True)
    return out/'bend_catchment_lookup.csv.gz'

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]
    out=root/'results/global_confinement/catchment_links'
    submission=Path('/Volumes/PhD/confinement/results/submission/v1')
    build_lookup(submission/'global_clusters_bend.gpkg',submission/'vector',out/'hydrobasins',out)
