

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 10:36:54 2024

@author: niekcollotdescury
"""
# import packages
from dataclasses import dataclass
import numpy as np
import shapely
import geopandas as gpd
import gc

# import partial packages
from shapely.geometry import LineString

# import custom modules
from .line_functions import create_angled_lines, adjust_confinement_line
from .extract_slope_along_raster_line import extract_slope_along_raster_line
from .calc_functions import extend_apex
from .dem import get_raster_vrt


MISSING_VALUE = 99999
RASTER_EXPANSION_FACTOR = 1.2


@dataclass
class OrthogonalGeometry:
    line_out: np.ndarray
    line_inn: np.ndarray
    left_right: np.ndarray


def remove_trailing_missing(arr, missing):
    arr  = np.asarray(arr)  # Ensure it's a NumPy array
    arr[arr == missing] = np.nan
    mask = np.isnan(arr)   # Identify missing values

    if not np.any(mask):   # No missing values
        return arr
    elif np.all(mask): # all missing valyes
        return np.nan
    else: # partial missing values
        last_non_missing_idx = np.where(~mask)[0][-1]  
        if np.any(mask[:last_non_missing_idx+1]): # missing not at the end
            return np.nan
        else: # only missing values at the end
            return arr[:last_non_missing_idx + 1]

def get_slope_values(slope_line, raster, slope_samples, demFillValeu):
                
    slope_line_angleP, slope_line_angleM = create_angled_lines(slope_line, 50)

    D, P     = extract_slope_along_raster_line(raster, 
                                                    slope_line, 
                                                    samples = slope_samples) 

    D_P, P_P = extract_slope_along_raster_line(raster, 
                                                slope_line_angleP, 
                                                samples = slope_samples) 

    D_M, P_M = extract_slope_along_raster_line(raster, 
                                                slope_line_angleM, 
                                                samples = slope_samples) 
    
    # Mean height profile line to mitigate dem measurement errors
    
    
    elevationValues = []
    for elev in [P, P_P, P_M]:

        elev = remove_trailing_missing(elev, demFillValeu)
        if isinstance(elev, np.ndarray):
            elevationValues.append(elev)

    if len(elevationValues) == 0:
        meanP, P, P_P, P_M,D, D_P, D_M = MISSING_VALUE, MISSING_VALUE, MISSING_VALUE, MISSING_VALUE, MISSING_VALUE, MISSING_VALUE, MISSING_VALUE
    else:
        min_length = min(len(arr) for arr in elevationValues)

        # Trim all arrays to the shortest length
        elevationValues = [arr[:min_length] for arr in elevationValues]
        D               = D[:min_length]
        meanP = np.mean(elevationValues, axis = 0)
    
    return (meanP, P, P_P, P_M, D, D_P, D_M)


def _empty_sampling_outputs():
    return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan


def _empty_orthogonal_outputs():
    return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan


def _format_profile_output(mean_profile, distances):
    if isinstance(mean_profile, int):
        elevation_values = MISSING_VALUE
    else:
        elevation_values = list(mean_profile)

    if isinstance(distances, int):
        distance_values = MISSING_VALUE
    else:
        distance_values = list(distances)

    return elevation_values, distance_values


def _open_reach_raster(df, bend_widths, reach_crs, dem_projection, cross_distance, dem_vrt, max_cross_distance):
    max_raster_size = max_cross_distance * RASTER_EXPANSION_FACTOR
    current_raster_size = min(np.max(bend_widths) * cross_distance * RASTER_EXPANSION_FACTOR, max_raster_size)
    raster = get_raster_vrt(dem_vrt, df, current_raster_size, reach_crs, dem_projection)
    if isinstance(raster, float):
        return np.nan, np.nan, current_raster_size, max_raster_size

    raster_border = shapely.box(*raster.rio.bounds()).exterior
    return raster, raster_border, current_raster_size, max_raster_size


def _expand_raster_to_cover_anchor(
    raster,
    raster_border,
    current_raster_size,
    max_raster_size,
    anchor_point,
    required_raster_size,
    df,
    reach_crs,
    dem_projection,
    dem_vrt,
):
    while anchor_point.distance(raster_border) < required_raster_size:
        if current_raster_size >= max_raster_size:
            break

        current_raster_size = min(current_raster_size * RASTER_EXPANSION_FACTOR, max_raster_size)
        raster.close()
        raster = get_raster_vrt(dem_vrt, df, current_raster_size, reach_crs, dem_projection)
        if isinstance(raster, float):
            return np.nan, np.nan, current_raster_size
        raster_border = shapely.box(*raster.rio.bounds()).exterior

    return raster, raster_border, current_raster_size


def _calculate_centerline_slope(raster, line, slope_samples):
    line_dist, line_elev = extract_slope_along_raster_line(raster, line, samples=slope_samples)

    line_dist = np.array(line_dist)
    line_elev[(line_elev < -300) | (line_elev == MISSING_VALUE)] = np.nan
    missing_elev = np.isnan(line_elev)
    missing_changes = np.diff(missing_elev).sum()

    line_dist = line_dist[missing_elev == False]
    line_elev = line_elev[missing_elev == False]

    if (missing_changes > 2) | (len(line_elev) == 0):
        return np.nan
    return np.polyfit(line_dist, line_elev, 1)[0]


def build_orthogonal_lines(
    line,
    apex_points,
    apexO_points,
    amplitudes,
    inf_lines,
    bend_widths,
    cross_distance,
):
    if (line.length <= 0) or isinstance(apex_points, float) or isinstance(inf_lines, float):
        return None

    loop_size = len(apex_points)
    line_out = np.empty(loop_size, dtype=object)
    line_inn = np.empty(loop_size, dtype=object)
    left_right = np.empty(loop_size)
    offset_line = line.offset_curve(30)

    for i in range(loop_size):
        apex_point = line.interpolate(line.project(apex_points[i]))
        apex_origin = apexO_points[i]
        bend_width = bend_widths[i]

        point_out, point_inn = extend_apex(
            apex_point,
            apex_origin,
            inf_lines[i],
            bend_width,
            amplitudes[i],
            cross_distance,
        )
        slope_line_out = adjust_confinement_line(LineString([apex_point, point_out]), line)
        slope_line_inn = adjust_confinement_line(LineString([apex_point, point_inn]), line)

        line_out[i] = slope_line_out
        line_inn[i] = slope_line_inn
        left_right[i] = 1 if slope_line_out.intersects(offset_line) else 0

    return OrthogonalGeometry(line_out=line_out, line_inn=line_inn, left_right=left_right)


def build_orthogonal_stage_frame(
    combined_reach_id,
    reach_crs: str,
    centerline,
    bend_lines,
    bend_widths,
    orthogonal_geometry: OrthogonalGeometry | None,
):
    if orthogonal_geometry is None:
        return None

    rows = [
        {
            "combined_reach_id": combined_reach_id,
            "bend_index": -1,
            "line_type": "centerline",
            "reach_crs": reach_crs,
            "bend_width": np.nan,
            "left_right": np.nan,
            "geometry": centerline,
        }
    ]

    for bend_index, bend_line in enumerate(bend_lines):
        bend_width = bend_widths[bend_index]
        left_right = orthogonal_geometry.left_right[bend_index]
        rows.extend(
            [
                {
                    "combined_reach_id": combined_reach_id,
                    "bend_index": bend_index,
                    "line_type": "bendline",
                    "reach_crs": reach_crs,
                    "bend_width": bend_width,
                    "left_right": left_right,
                    "geometry": bend_line,
                },
                {
                    "combined_reach_id": combined_reach_id,
                    "bend_index": bend_index,
                    "line_type": "outer",
                    "reach_crs": reach_crs,
                    "bend_width": bend_width,
                    "left_right": left_right,
                    "geometry": orthogonal_geometry.line_out[bend_index],
                },
                {
                    "combined_reach_id": combined_reach_id,
                    "bend_index": bend_index,
                    "line_type": "inner",
                    "reach_crs": reach_crs,
                    "bend_width": bend_width,
                    "left_right": left_right,
                    "geometry": orthogonal_geometry.line_inn[bend_index],
                },
            ]
        )

    return gpd.GeoDataFrame(rows, geometry="geometry", crs=reach_crs)


def sample_orthogonal_profiles(
    line,
    df: 'gpd.GeoDataSeries',
    bend_lines,
    bend_widths,
    orthogonal_geometry: OrthogonalGeometry | None,
    reach_crs: str,
    dem_projection,
    cross_distance: int,
    demFillValeu,
    demVRT,
    slope_samples=100,
    maxCrossDistance=30000,
):
    if orthogonal_geometry is None:
        return _empty_sampling_outputs()

    raster, raster_border, current_raster_size, max_raster_size = _open_reach_raster(
        df,
        bend_widths,
        reach_crs,
        dem_projection,
        cross_distance,
        demVRT,
        maxCrossDistance,
    )
    if isinstance(raster, float):
        return _empty_sampling_outputs()

    loop_size = len(orthogonal_geometry.line_out)
    elev_out, elev_inn = np.empty(loop_size, dtype=object), np.empty(loop_size, dtype=object)
    dist_out, dist_inn = np.empty(loop_size, dtype=object), np.empty(loop_size, dtype=object)
    bend_height = np.empty(loop_size)

    try:
        for i in range(loop_size):
            bend_width = bend_widths[i]
            anchor_point = shapely.Point(orthogonal_geometry.line_out[i].coords[0])
            required_raster_size = bend_width * cross_distance if bend_width * cross_distance < maxCrossDistance else maxCrossDistance

            raster, raster_border, current_raster_size = _expand_raster_to_cover_anchor(
                raster,
                raster_border,
                current_raster_size,
                max_raster_size,
                anchor_point,
                required_raster_size,
                df,
                reach_crs,
                dem_projection,
                demVRT,
            )
            if isinstance(raster, float):
                return _empty_sampling_outputs()

            _, bend_height_profile = extract_slope_along_raster_line(raster, bend_lines[i], samples=slope_samples)

            mean_out, _, _, _, distances_out, _, _ = get_slope_values(
                orthogonal_geometry.line_out[i],
                raster,
                slope_samples,
                demFillValeu,
            )
            mean_inn, _, _, _, distances_inn, _, _ = get_slope_values(
                orthogonal_geometry.line_inn[i],
                raster,
                slope_samples,
                demFillValeu,
            )

            elev_out[i], dist_out[i] = _format_profile_output(mean_out, distances_out)
            elev_inn[i], dist_inn[i] = _format_profile_output(mean_inn, distances_inn)

            if np.isnan(np.nanmean(bend_height_profile)):
                bend_height[i] = MISSING_VALUE
            else:
                bend_height[i] = np.nanmean(bend_height_profile)

        line_slope = _calculate_centerline_slope(raster, line, slope_samples)
    finally:
        if hasattr(raster, "close"):
            raster.close()
        del raster, raster_border
        gc.collect()

    return elev_out, elev_inn, dist_out, dist_inn, line_slope, bend_height

def get_orthogonals(line,df:'gpd.GeoDataSeries',
                    apex_points, apexO_points, amplitudes, infLines,
                    bendLines, bendWidths,
                    reachCRS:str, DEMprojection,
                    cross_distance:int,
                    demFillValeu, demVRT,
                    directory:str, slope_samples = 100,
                    maxCrossDistance = 30000):
    """
    Function that creates lines at approximately right angles with the river 
    and determines the slope of these lines.\n
    input:
    - line: centerline LineString
    - df: Reach Dataframe (DataSeries)
    - Raster: DEM raster clipped to river segment extent
    - apex_points: list of apex points in local projeciton
    - apexO_points: list of apex origin points on the inflection line
    - inf_points: list of inflection points
    - bendLines: list of bendLines
    - bendWidths: list of bendWidths
    - reachCRS: local utm projection
    - DEMprojection: projection of DEM
    - width_input: reach width --> check for raster size
    - cross_distance: the factor that half the width is multiplied with to determine the length of 
                      the line created at right angels with the river. 
    - demFillValue: fill value for the selected DEM
    - demVRT: virtual dataset for DEM
    - directory: folder directory
    - slope_samples: Number of samples extracted in slope function\n
    
    return:
    - Elevation profile outerbend
    - Elevation profile innerbend
    - distance outbend
    - distance innerbend 
    - linestring of outerbend
    - linestring of innerbend
    - binary for matching left with outer bend. 1 is left, 0 is right
    - centerlineslope\n

    If confining slope line intersects with the centerline the line is cutoff at the centerline. 
    This can create differences in length between the inner and outer bend
    """
    orthogonal_geometry = build_orthogonal_lines(
        line,
        apex_points,
        apexO_points,
        amplitudes,
        infLines,
        bendWidths,
        cross_distance,
    )
    if orthogonal_geometry is None:
        return _empty_orthogonal_outputs()

    elev_out, elev_inn, dist_out, dist_inn, line_slope, bend_height = sample_orthogonal_profiles(
        line,
        df,
        bendLines,
        bendWidths,
        orthogonal_geometry,
        reachCRS,
        DEMprojection,
        cross_distance,
        demFillValeu,
        demVRT,
        slope_samples=slope_samples,
        maxCrossDistance=maxCrossDistance,
    )

    return (
        elev_out,
        elev_inn,
        dist_out,
        dist_inn,
        orthogonal_geometry.line_out,
        orthogonal_geometry.line_inn,
        orthogonal_geometry.left_right,
        line_slope,
        bend_height,
    )
