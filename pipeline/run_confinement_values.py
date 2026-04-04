from .bend_io import read_bend_table
from .support import confinement_factor, str_to_list, str_to_list_comb
from .calc_functions import confinement_values
from .paths import load_project_paths, format_factor_token

import numpy as np
import os
import pandas as pd
from pathlib import Path

from tqdm import tqdm 


def _split_file_stem(file_stem):
    parts = str(file_stem).split("_")
    if len(parts) >= 2:
        return parts[0], parts[1]
    return None, None


def _list_length(value):
    if isinstance(value, list):
        return len(value)
    return 0


def _row_bend_count(row, candidate_cols):
    for col in candidate_cols:
        length = _list_length(row.get(col))
        if length > 0:
            return length
    return 0


def _value_for_bend(value, bend_index):
    if isinstance(value, list):
        if bend_index < len(value):
            return value[bend_index]
        return np.nan
    return value


def prepare_bend_level_dataframe(df, file_stem=None):
    dfInc = df[(df['include_flag'] == '0') & (df['calculated'] == '0')].copy()
    if dfInc.empty:
        df_empty = dfInc.copy()
        if 'bend_index' not in df_empty.columns:
            df_empty['bend_index'] = pd.Series(dtype='int64')
        if file_stem is not None:
            file_cont, file_id = _split_file_stem(file_stem)
            df_empty['file_stem'] = file_stem
            df_empty['file_cont'] = file_cont
            df_empty['file_id'] = file_id
            df_empty['bend_id'] = pd.Series(dtype='object')
        return df_empty

    dropCols = ['x', 'y', 'swot_obs', 'edit_flag', 'trib_flag', 'wse_var', 'width_var']
    dfInc = dfInc.drop(columns=[col for col in dropCols if col in dfInc.columns])

    listStrCols = ['LROrthog', 'bendMaxWidths', 'bendWidths', 'ang', 'bendSin', 'apex', 'bendDistOut', 'bendLen', 'bendHeight']
    listGeomCols = ['apexP', 'lineInn', 'lineOut', 'bendLines']
    listNestCols = ['distOut', 'distInn', 'elevInn', 'elevOut']
    listCols = listStrCols + listGeomCols + listNestCols

    dfReach = str_to_list_comb(dfInc, listCols, listNestCols)

    if 'geometry' in dfReach.columns:
        dfReach['geometry'] = dfReach['geometry'].apply(
            lambda geom: geom.wkt if hasattr(geom, 'wkt') else geom
        )

    bend_count_cols = ['bendWidths', 'bendLines', 'distOut', 'elevOut', 'apex', 'lineInn', 'lineOut', 'LROrthog']
    bend_rows = []
    for _, row in dfReach.iterrows():
        bend_count = _row_bend_count(row, bend_count_cols)
        if bend_count == 0:
            continue

        row_dict = row.to_dict()
        for bend_index in range(bend_count):
            bend_row = row_dict.copy()
            for col in listCols:
                bend_row[col] = _value_for_bend(row_dict.get(col), bend_index)
            bend_row['bend_index'] = bend_index
            bend_rows.append(bend_row)

    dfE = pd.DataFrame(bend_rows)
    if dfE.empty:
        dfE = pd.DataFrame(columns=list(dfReach.columns) + ['bend_index'])

    if file_stem is not None:
        file_cont, file_id = _split_file_stem(file_stem)
        dfE['file_stem'] = file_stem
        dfE['file_cont'] = file_cont
        dfE['file_id'] = file_id
        dfE['bend_id'] = (
            dfE['file_stem'].astype(str)
            + '_'
            + dfE['combined_reach_id'].astype(int).astype(str)
            + '_'
            + dfE['bend_index'].astype(int).astype(str)
        )

    return dfE


def _format_height_factor(hf):
    if isinstance(hf, float) and hf.is_integer():
        hf = int(hf)
    hf_save = f"{hf}"
    if hf < 10:
        hf_save = f"0{hf}"
    return hf_save


def _resolve_step5_paths(config_path=None, conf_factor=50):
    project_paths = load_project_paths(config_path)
    project_paths.ensure_step5_dirs()

    return {
        "project_paths": project_paths,
        "results_root": project_paths.results_root,
        "single_values_dir": project_paths.single_values_dir,
        "reach_averaged_dir": project_paths.reach_averaged_dir,
        "confinement_factor_file": project_paths.confinement_factor_file_for(conf_factor),
    }


def _resolve_step6_paths(config_path=None):
    project_paths = load_project_paths(config_path)
    project_paths.ensure_step5_dirs()

    return {
        "project_paths": project_paths,
        "results_root": project_paths.results_root,
        "single_values_dir": project_paths.single_values_dir,
        "reach_averaged_dir": project_paths.reach_averaged_dir,
    }


def ER_slope_margin_values(dfInc,demVRT, cf = [50,10], heightFactor = 2,
                           config_path = None):
    from scipy.spatial import KDTree
    import geopandas as gpd
    import shapely

    print(f'run_confinement_values - calc_confinement_values - ER_slope_margin_values: Start')
    
    # Find confinement factor value
    step5_paths = _resolve_step5_paths(config_path, cf[0])
    confinement_factor_file = step5_paths["confinement_factor_file"]
    if confinement_factor_file.exists() is False:
        raise FileNotFoundError(
            f"Confinement factor table not found at {confinement_factor_file}. "
            "Run 'python -m pipeline.build_step4_confinement_factor' first."
        )
    dfCF = pd.read_csv(confinement_factor_file)
    # Using KDTree for efficient nearest neighbor search
    tree = KDTree(dfCF[['bendWidths']].values)
    _, idx = tree.query(dfInc[['bendWidths']].values)
    # Getting the closest values
    dfInc['conFactor'] = dfCF.iloc[idx]['conFactor'].values

    
    # change to combined reach id!!!!!!!! so that local crs is correct
    # reachIDS = dfInc['reach_id'].values
    reachIDS = dfInc['combined_reach_id'].astype(int).unique()
    crGeoms  = np.empty(len(reachIDS), dtype = object)
    # for i, rid in tqdm(enumerate(reachIDS), total = len(reachIDS)):
    for i, rid in enumerate(reachIDS):

        dfReach     = dfInc[dfInc['combined_reach_id'] == rid].copy()
        centerline_vals = dfReach.get('centerlineWkt')
        if centerline_vals is not None:
            centerline_vals = centerline_vals.dropna()
        if centerline_vals is not None and centerline_vals.empty is False:
            crGeoms[i] = shapely.from_wkt(centerline_vals.iloc[0])
        else:
            print(f"Missing centerlineWkt for combined_reach_id {rid}; saving empty reach-averaged geometry.")
            crGeoms[i] = shapely.geometry.LineString()


        for i, r in dfReach.iterrows():
            dfBend = gpd.GeoDataFrame({'geometry':shapely.from_wkt(r['bendLines'])}
                                      , index = [0], crs = r['combined_reach_crs'])
            dfBend   = dfBend.to_crs('EPSG:4326')
            bendGeom = dfBend['geometry'].iloc[0]

            (vwo, vwi, cph, cmh, so, si, 
                ERo, ERi) = confinement_values(r['elevOut'], r['elevInn'],
                                                r['distOut'], r['distInn'], 
                                                r['bendWidths'],r['bendMaxWidths'],
                                                r['conFactor'] * heightFactor)
            
            dfInc.loc[i, ['valley_width_out', 'valley_width_inn', 
                    'cp_height', 'cm_height', 
                    'slope_out', 'slope_inn',
                    'ER_out', 'ER_inn', 'bendGeom'
                    ]] = (vwo, vwi, cph, cmh, so, si, ERo, ERi, bendGeom)



    dfInc['ER_left']  = np.where(dfInc['LROrthog'] == 1, dfInc['ER_out'], dfInc['ER_inn'])
    dfInc['ER_right'] = np.where(dfInc['LROrthog'] == 1, dfInc['ER_inn'], dfInc['ER_out'])

    dfInc['slope_left']  = np.where(dfInc['LROrthog'] == 1, dfInc['slope_out'], dfInc['slope_inn'])
    dfInc['slope_right'] = np.where(dfInc['LROrthog'] == 1, dfInc['slope_inn'], dfInc['slope_out'])

    return dfInc, crGeoms

def calc_confinement_values(df,fileName, returnDataframe, open_seperate = False,
                            crossFactor = 50, hf = 2, config_path = None):
    from osgeo import gdal
    import geopandas as gpd

    print(f'run_confinement_values - calc_confinement_values: {fileName}')
    gdal.UseExceptions()
    step5_paths = _resolve_step5_paths(config_path, crossFactor)
    project_paths = step5_paths["project_paths"]
    vrt_file = str(project_paths.fabdem_vrt)
    demVRT   = gdal.Open(vrt_file)
    if demVRT is None:
        raise FileNotFoundError(
            f"FABDEM VRT not found at {vrt_file}. Run 'python -m pipeline.build_fabdem_index' first."
        )

    hfSave = _format_height_factor(hf)

    if open_seperate == True:
        listCols = ['distOut', 'distInn', 'elevInn', 'elevOut']
        df = str_to_list(df, listCols, listCols)

    dfE, crGeoms = ER_slope_margin_values(
        df,
        demVRT,
        [50,10],
        hf,
        config_path=config_path,
    )

    ########################################
    # Save reach averaged gpkg
    ########################################
    dfEG = dfE.groupby('combined_reach_id').agg({
    'ER_out':'mean'    ,'ER_inn':'mean',
    'ER_left':'mean'   ,'ER_right':'mean',
    'slope_out':'mean' ,'slope_inn': 'mean',
    'slope_left':'mean','slope_right': 'mean',
    'cp_height':'mean',
    'catchment_position':'first'}).copy()


    gdfEG = gpd.GeoDataFrame(dfEG, geometry = crGeoms, crs = 'EPSG:4326')
    reach_averaged_file = step5_paths["reach_averaged_dir"] / f"{fileName}_{hfSave}.gpkg"
    gdfEG.to_file(reach_averaged_file, driver = 'GPKG')
    
    ########################################
    # Save reduced bend gpkg
    ########################################
    saveCols = ['reach_id', 'combined_reach_id',
                'apex', 'bendWidths', 'bendMaxWidths',
                'ang', 'sin', 'bendSin','cp_height',
                'ER_out', 'ER_inn','ER_left', 'ER_right',
                'slope_out', 'slope_inn', 'slope_left', 'slope_right',
                'bendDistOut', 'bendLen']
    dfE_bend = gpd.GeoDataFrame(dfE[saveCols].copy(), geometry = dfE['bendGeom'], crs= 'EPSG:4326')
    bend_gpkg_file = step5_paths["single_values_dir"] / f"{fileName}_{hfSave}_conf.gpkg"
    dfE_bend.to_file(bend_gpkg_file, driver = 'GPKG')


    dfE = dfE.drop(['elevOut', 'elevInn', 'distOut', 
                    'distInn', 'geometry', 'bendGeom'], axis = 1)
    dfE['apex'] = dfE['apex'].astype('float64')
    # save nc File
    dfIncX = dfE.to_xarray()
    ncFile = step5_paths["single_values_dir"] / f"{fileName}_{hfSave}_conf.nc"
    if ncFile.exists():
        ncFile.unlink()
    
    dfIncX.to_netcdf(ncFile)

    if returnDataframe == True:
        return dfE

def concat_nc_conf_files(cross = 50, hf = 2, config_path = None):
    import xarray as xr

    step6_paths = _resolve_step6_paths(config_path)
    cross_token = format_factor_token(cross)
    hf_token = _format_height_factor(hf)
    dsList = []
    
    for c in ['af', 'as', 'sa', 'na', 'oc', 'eu']:
        print(c)
        files = sorted(step6_paths["single_values_dir"].glob(f"{c}_??_{cross_token}_{hf_token}_conf.nc"))
        for f in tqdm(files):
            dsTemp = xr.open_dataset(f)
            dsTemp['file'] = ('index', [c] * dsTemp.sizes['index'])
            dsTemp = dsTemp.drop_vars(['infP', 'bendLines', 'apexP', 'lineInn', 'lineOut'])

            parts = Path(f).stem.split("_")
            dsTemp['file_cont'] = parts[0]
            dsTemp['file_num']  = parts[1]

            dsList.append(dsTemp)

    if len(dsList) == 0:
        raise FileNotFoundError(
            f"No Step 5 NetCDF outputs were found in {step6_paths['single_values_dir']} "
            f"for cross {cross_token} and height factor {hf_token}."
        )

    ds = xr.concat(dsList, dim='index')
    output_file = step6_paths["single_values_dir"] / f"global_{cross_token}_{hf_token}_conf.nc"
    if output_file.exists():
        output_file.unlink()
    ds.to_netcdf(output_file)
    return output_file

def concat_reachAveraged(cross = 50, hf = 2, config_path = None):
    import geopandas as gpd

    step6_paths = _resolve_step6_paths(config_path)
    cross_token = format_factor_token(cross)
    hf_token = _format_height_factor(hf)
    output_files = []
    for c in ['af', 'sa', 'as', 'na', 'eu', 'oc']:
        print(c, cross_token, hf_token)

        files = sorted(step6_paths["reach_averaged_dir"].glob(f"{c}_??_{cross_token}_{hf_token}.gpkg"))
        if len(files) == 0:
            continue

        for i, f in enumerate(files):
            parts = Path(f).stem.split("_")
            contN = parts[0]
            contI = parts[1]
            dfT = gpd.read_file(f)
            dfT['file_cont'] = contN
            dfT['file_id']   = contI
            if i == 0:
                dfW = dfT
            else:
                dfW = pd.concat([dfW, dfT])
        output_file = step6_paths["reach_averaged_dir"] / f"{c}_{cross_token}_{hf_token}_reachAveraged_conf.gpkg"
        if output_file.exists():
            output_file.unlink()
        dfW.to_file(output_file, driver = 'GPKG')
        output_files.append(output_file)

    if len(output_files) == 0:
        raise FileNotFoundError(
            f"No Step 5 reach-averaged outputs were found in {step6_paths['reach_averaged_dir']} "
            f"for cross {cross_token} and height factor {hf_token}."
        )
    return output_files


def create_apex_val_dataframe(file, output_path = None, returnDataframe = False, config_path = None):
    paths = load_project_paths(config_path)
    paths.ensure_step3_dirs()
    input_file = Path(file).expanduser()
    if input_file.is_absolute() is False:
        input_file = (paths.repo_root / input_file).resolve()

    if output_path is None:
        output_file = paths.single_values_dir / input_file.name
    else:
        output_file = Path(output_path).expanduser()
        if output_file.is_absolute() is False:
            output_file = (paths.repo_root / output_file).resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if input_file.suffix.lower() == ".parquet":
        dfE = read_bend_table(input_file)
    else:
        df = pd.read_csv(input_file, dtype = {'include_flag': str, 'calculated': str})
        dfE = prepare_bend_level_dataframe(df, file_stem=input_file.stem)
    
    print('Start saving:', output_file.name)
    dfE.to_csv(output_file, index=False)

    # calc_confinement_values(dfE, file[-12:-4], directory, returnDataframe)
    if returnDataframe == True:
        return dfE
