# import packages
import warnings
import shapely

#import partial packages
from shapely.geometry import LineString, Point

# import packages for function typecasting
import geopandas as gpd
from typing import List, Tuple

def connect_multilinestring_rows(rows:gpd.GeoDataFrame) -> LineString:
    """
    Function for connecting the coordinates of a multilinestring in the correct order.\n
    input:
    - rows: pandas dataseries with the rows that have to be connected\n
    output:
    - newGeom: LineString of connected coordinates
    """
    rows = rows.sort_values('reach_order')
    rid  = list(rows.reach_id.values)

    newGeom, ID, conID    = [],  rid[0], rows.iloc[0]['dn_connected_reach']
    distEval = 'd1<d2'
    for i in range(len(rid)):
        if i > 0:
            conID = ID
            ID    = newID
            distEval = 'd1>d2'
        row   = rows[rows.reach_id == ID].iloc[0]
        newID = row['dn_connected_reach']

        rowGeom = rows[rows.reach_id == ID].iloc[0].geometry
        conGeom = rows[rows.reach_id == conID].iloc[0].geometry
        d1 = Point(rowGeom.coords[0]).distance(conGeom)
        d2 = Point(rowGeom.coords[-1]).distance(conGeom)
        if  eval(distEval):
            newGeom.extend(rowGeom.coords[::-1])
        else:
            newGeom.extend(rowGeom.coords)
    newGeom = LineString(newGeom)
    return newGeom

def check_line_start(reachLine, dfReach, df, localCRS):
    """
    Check the start of line coordinates.
    Make sure the start is connected to upstream reach.\n
    input:
    - dfReach: reaches in combined reach id
    - df: complete dataframe
    - localCRS: dfReach projection\n
    output:
    - reachLine: updates lineString
    """
    dfReach = dfReach.copy()
    df      = df.copy()


    def flip(line, conID, df):
        conRow = df[df['reach_id'] == conID]
        conRow = conRow.to_crs(localCRS)
        conLine = conRow.iloc[0].geometry

        distStart = Point(line.coords[0]).distance(conLine)
        distEnd   = Point(line.coords[-1]).distance(conLine)
        if distStart > distEnd:
            line = LineString(list(line.coords)[::-1])

        
        return line
    
    rchUP = dfReach.loc[dfReach['reach_order'] == 0, 'rch_id_up_orig'].iloc[0]
    rchDN = dfReach.loc[dfReach['reach_order'] == dfReach['reach_order'].max(), 'rch_id_dn_orig'].iloc[0]

    if isinstance(rchUP, list):
        reachLine = flip(reachLine, rchUP[0], df)
    elif isinstance(rchDN, list):
        reachLine = flip(reachLine, rchDN[0], df)
    else:
        print(f'Connect_geometries - merge_centerlines - check_line_start: LoneReach ({dfReach.combined_reach_id.unique()})')
    return reachLine

def merge_centerlines(dfLines:gpd.GeoDataFrame, df:gpd.GeoDataFrame, 
                      localCRS:str, adjustStartSide = True)-> Tuple[LineString, float, bool]:
    """
    Merge seperate adjoining linestrings into single linestring.\n
    input: 
    - dfLines: geopandas GeoDataFrame consisting of the rows to be connected
    - df: total dataframe (including skipped reaches
    - localCRS: reach crs
    - adjustStartSide: boolean to change the side of the linestring to upstream\n
    output:
    - line: new linestring
    - line.length: new line length
    - wrongDiff: boolean value if the difference between the sum of the seperate linestrings has an absolute difference with the connected linestrings of more than 200m
    """
    lines   = dfLines.geometry.values

    if dfLines.shape[0] > 1:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            multiLines = shapely.unary_union(lines)
            singleLine = shapely.line_merge(multiLines)

        if singleLine.geom_type == 'LineString':
            line = singleLine
            A = 1
        elif singleLine.geom_type == 'MultiLineString':
            # remove self intersecting lines
            A =2
            newLine = []
            for l in singleLine.geoms:          
                if len(shapely.ops.polygonize(l)) ==0:
                    newLine.append(l)
            line = shapely.line_merge(shapely.unary_union(newLine))
            
            # backup connection of linestrings
            if line.geom_type != 'LineString':
                A = 3
                line = connect_multilinestring_rows(dfLines)
        else: 
            print(f'merge_centerlines: wrong geom type for lines input: {dfLines['combined_reach_id'].unique()}')

    else:
        line = lines[0]
    
    lineLen = 0
    for l in lines:
        lineLen += l.length

    lineDiff = abs(line.length - lineLen)
    wrongDiff = False
    if lineDiff > 200:
        wrongDiff = True
    if adjustStartSide == True:
        line = check_line_start(line, dfLines, df, localCRS)

    return line,line.length, wrongDiff