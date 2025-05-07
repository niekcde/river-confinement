# Import packages
import numpy as np
import shapely
import matplotlib.pyplot as plt
import geopandas as gpd

# Import modules from package
from scipy.spatial import cKDTree
from shapely.geometry import LineString, MultiLineString, Polygon, Point

#import custom modules
from line_functions import getExtrapoledLine
from extract_slope_along_raster_line import extract_slope_along_raster_line
from dem import get_raster_vrt

def get_raster_coordinates(raster):
    x = raster.x.values
    y = raster.y.values
    xx, yy = np.meshgrid(x, y)
    coords = np.column_stack([xx.ravel(), yy.ravel()])
    return coords

# Find the closest point with a specific elevation relative to the line
def find_raster_points_threshold(dem, target_elevation):
    # Get DEM coordinates
    dem_coords = get_raster_coordinates(dem)
    
    # Get elevation data
    dem_elevations = dem.values.ravel()
    # Filter out points that are not at the target elevation
    valid_points = np.where(dem_elevations > target_elevation)[0]

    if len(valid_points) == 0:
        raise ValueError(f"No points with elevation {target_elevation} found.")
    
    valid_coords = dem_coords[valid_points]
    
    # Use cKDTree for fast nearest neighbor search
    tree = cKDTree(valid_coords)
    return tree, valid_coords

def get_raster_side(line, demVRT, bufferSize, localCRS):
    if isinstance(line, str):
        line = shapely.from_wkt(line)
    # BSize      = row['bendMaxWidths'] * (distanceFactor + 0.1)
    # BSize      = row['bendMaxWidths'] * (distanceFactor + 0.1)

    buffer     = line.buffer(bufferSize)
    if buffer.geom_type == 'MultiPolygon':
        BL = 0
        for b, B in enumerate(buffer.geoms):
            if B.length > BL:
                BL = B.length
                Bind = b
        buffer = buffer.geoms[Bind]
    bufferLine = LineString(buffer.exterior.coords)

    p1, p2 = line.interpolate(50), line.interpolate(150)
    p3, p4 = line.interpolate(-50), line.interpolate(-150)

    start    = getExtrapoledLine(p2.coords[0], p1.coords[0], 
                                (bufferSize + 150) * 1.5, 1, True, 'Distance')
    end      = getExtrapoledLine(p4.coords[0], p3.coords[0],
                                (bufferSize+ 150) * 1.5, 1, True, 'Distance')

    if start.intersects(bufferLine) == False:
        bufferPointS = shapely.ops.nearest_points(Point(start.coords[-1]), bufferLine)[1]
        start        = LineString([p1, bufferPointS])
    if end.intersects(bufferLine) == False:
        bufferPointS = shapely.ops.nearest_points(Point(end.coords[-1]), bufferLine)[1]
        start        = LineString([p3, bufferPointS])

    startEnd = MultiLineString([start, end])
    # Split buffer in two seperate sections
    bufferLineSplit = shapely.ops.split(bufferLine, startEnd).geoms 

    ##### add assertion statement about number of splits!
    comb =[]
    for i in range(len(bufferLineSplit)):
        bufferSection  = bufferLineSplit[i]
        d1 = bufferSection.distance(start) 
        d2 = bufferSection.distance(end)

        if (d1 > 1) | (d2 > 1):
            comb.append(bufferSection)
        else:
            bufferOutside1 = bufferSection
    bufferOutside2 = shapely.line_merge(MultiLineString(comb))

    def create_poly_from_lines(L1, L2):

        P1 = Polygon([*list(L1.coords), *list(L2.coords)[::-1]])
        P2 = Polygon([*list(L1.coords), *list(L2.coords)])

        if P1.exterior.length > P2.exterior.length:
            P = P2
        else:
            P = P1
        return P

    bufferSide1 = create_poly_from_lines(bufferOutside1, line)
    bufferSide2 = create_poly_from_lines(bufferOutside2, line)


    ###################################
    # Open, reproject and clip dem raster
    ###################################

    dfLine = gpd.GeoDataFrame({'geometry':line}, index =[0], crs = localCRS)
    raster = get_raster_vrt(demVRT,dfLine ,bufferSize,localCRS, 'EPSG:4326')
    raster = raster.rio.write_crs(localCRS)

    # clip raster to buffer sides
    rasterSide1 = raster.rio.clip([bufferSide1])
    rasterSide2 = raster.rio.clip([bufferSide2])
    return raster, rasterSide1, rasterSide2

def confinement_margin(row,raster, raster1, raster2, distanceFactor = 2, heightFactor = 2, stepsize = 50):
    """
    input:
    - reach: SWORD Row, GeoDataFrame
    - df_node: SWORD Nodes, GeoDataFrame
    - directory: code directory
    - projection: global projection (default: global mercator)
    - distanceFactor: Maximum width multiplier from the centerline for confinement
    - heightFactor: Height above river centerline
    - stepsize: confining distance step (50m)
    - plot: boolean for plotting (False)\n

    return:
    - Percentage of centerline that is confinement. Total and per side
    """

    ###################################
    # Create buffer around centerline and extend centerline to intersect with buffer
    #   intersection used to create left and right buffer sections
    ###################################
    line = row['bendLines']
    if isinstance(line, str):
        line = shapely.from_wkt(line)

    ###################################
    # Calculate distance from centerline to confining margin
        # Double confinement counts single
    ###################################
    def conf_distance(raster1, raster2, width, threshold, rivHeight, line, confStepSize, distanceFactor):
        def find_cond_distance(line, confStepSize, tree, coords, distanceFactor):
            distance = 0
            p, lineIndex, TP = [], [], []
            step = confStepSize / line.length
            for i, I in enumerate(np.arange(step,1 + step,step)):
                targetPoint = line.interpolate(I, normalized = True)
                TQ = tree.query(targetPoint.coords[0])
                TP.append(targetPoint)
                widthThreshold = width*distanceFactor
                cfs = confStepSize
                if I > 1:
                    cfs = line.length - (i) * confStepSize
                    
                if (TQ[0] < widthThreshold):
                    distance += cfs
                    p.append(coords[TQ[1]])
                    lineIndex.append(i)

            return distance,lineIndex, p, TP, cfs

        ind1, ind2, p1, p2, TP = [],[], [], [], []
        d1, d2 = 0, 0
        if threshold < raster1.max():
            tree1, coords1        = find_raster_points_threshold(raster1, threshold)
            d1,ind1,  p1, tp, cfs = find_cond_distance(line, confStepSize, tree1, coords1,
                                                        distanceFactor)
            TP.extend(tp)

        if threshold < raster2.max():
            tree2, coords2       = find_raster_points_threshold(raster2, threshold)
            d2,ind2,  p2, _, cfs = find_cond_distance(line, confStepSize, tree2, coords2, distanceFactor)

        double_confinement  = list(set(ind1) & set(ind2))
        distance = (d1 + d2) - (len(double_confinement) * confStepSize)
        if (len(ind1) > 0) & (len(double_confinement) > 0):
            if double_confinement[-1] == ind1[-1]:
                distance = (d1 + d2) - (((len(double_confinement)-1) * confStepSize) + cfs)
        
        return distance, d1,d2, p1, p2, TP

    ###################################
    # Run code for each Bend
    ###################################
    line_values = extract_slope_along_raster_line(raster, line,400)
    lineHeight = np.median(line_values)


    threshold = (row['conFactor'] * heightFactor * row['bendWidths']) + lineHeight
    D, D1, D2, p1, p2,tp = conf_distance(raster1, raster2, row['bendWidths'], 
                                            threshold, lineHeight, line, stepsize, distanceFactor)


    # f, ax = plt.subplots()
    # ax.plot(*line.xy)
    # rasterValues = raster[:,:]
    # rasterValues = rasterValues.values.ravel()
    # vmin = rasterValues[rasterValues != -9999].min()
    # for p in p1:
    #     ax.scatter(*p, zorder = 1000, s = 20, c = 'red', marker = '^')

    # for p in p2:
    #     ax.scatter(*p, zorder = 1000, s = 20, c = 'red', marker = '^')

    # raster.plot(ax= ax, vmin = vmin)
    # plt.show()


    # adjust to confinement per bend!!!!!!!
    return ((D / line.length), (D1 / line.length), (D2 / line.length), lineHeight,
            p1, p2)