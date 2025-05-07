# Import packages
import numpy as np
import shapely
import matplotlib.pyplot as plt

# Import modules from package
from scipy.spatial import cKDTree
from shapely.geometry import LineString, MultiLineString, Polygon

#import custom modules
from support import get_ids
from line_functions import getExtrapoledLine
from extract_slope_along_raster_line import extract_slope_along_raster_line
from dem import find_dem

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


def confinement_margin(df_reach, df_node, directory:str,projection = 'EPSG:3857', distanceFactor = 2, heightFactor = 2, stepsize = 50, plot = False):
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
    # Filter dataframes and change crs to local crs 
    ###################################
    df_reach_in = df_reach.copy()
    df_reach    = df_reach.to_crs(df_reach.iloc[0].localCRS)
    reach       = df_reach.iloc[0]

    reachIds   = get_ids(reach)
    reachNodes = df_node[df_node.reach_id.isin(reachIds)]
    reachNodes = reachNodes.reset_index()
    reachNodes = reachNodes.to_crs(reach.localCRS)

    # sort the nodes along the centerline. So that index 0 is the start and -1 is the end of the centerline
    nodePosOnLine = [reach.geometry.project(n.geometry) for i, n in reachNodes.iterrows()]
    nodeIndex    = np.argsort(nodePosOnLine)
    reachNodes = reachNodes.iloc[nodeIndex]



    ###################################
    # Create buffer around centerline and extend centerline to intersect with buffer
    #   intersection used to create left and right buffer sections
    ###################################
    
    BSize      = reachNodes.max_width.max() * 5 # buffer size
    line       = reach.geometry
    buffer     = line.buffer(BSize)
    if buffer.geom_type == 'MultiPolygon':
        BL = 0
        for b, B in enumerate(buffer.geoms):
            if B.length > BL:
                BL = B.length
                Bind = b
        buffer = buffer.geoms[Bind]


    bufferLine = LineString(buffer.exterior.coords)

    node1, node2 = reachNodes.iloc[0 ].geometry, reachNodes.iloc[1 ].geometry, 
    node3, node4 = reachNodes.iloc[-1].geometry, reachNodes.iloc[-2].geometry
    start    = getExtrapoledLine(node2.coords[0], node1.coords[0], 
                                 (BSize + node1.distance(node2)) * 1.5, 1, True, 'Distance')
    end      = getExtrapoledLine(node4.coords[0], node3.coords[0],
                                  (BSize+ node3.distance(node4)) * 1.5, 1, True, 'Distance')
    
    if start.intersects(bufferLine) == False:

        bufferPointS = shapely.ops.nearest_points(node1, bufferLine)[1]
        start       = getExtrapoledLine(node1.coords[0], bufferPointS.coords[0],
                                node1.distance(bufferPointS) + 200, 1, True, 'Distance')

    if end.intersects(bufferLine) == False:

        bufferPointE = shapely.ops.nearest_points(node3, bufferLine)[1]
        end         = getExtrapoledLine(node3.coords[0], bufferPointE.coords[0],
                                node3.distance(bufferPointE) + 200, 1, True, 'Distance')


    startEnd = MultiLineString([start, end])
    # Split buffer in two seperate sections
    bufferLineSplit = shapely.ops.split(bufferLine, startEnd).geoms 

    ##### add assertion statement about number of splits!
    comb =[]
    for i in range(len(bufferLineSplit)):
        e  = bufferLineSplit[i]
        d1 = e.distance(start) 
        d2 = e.distance(end)

        if (d1 > 1) | (d2 > 1):
            comb.append(e)
        else:
            p1 = e
    p2 = shapely.line_merge(MultiLineString(comb))

    def create_poly_from_lines(L1, L2):

        P1 = Polygon([*list(L1.coords), *list(L2.coords)[::-1]])
        P2 = Polygon([*list(L1.coords), *list(L2.coords)])

        if P1.exterior.length > P2.exterior.length:
            P = P2
        else:
            P = P1
        return P

    bufferSide1 = create_poly_from_lines(p1, line)
    bufferSide2 = create_poly_from_lines(p2, line)


    ###################################
    # Open, reproject and clip dem raster
    ###################################
    demWidth = BSize * 2
    raster = find_dem(df_reach_in, directory, projection,
                        demWidth, False)  
    raster     = raster.rio.reproject(reach.localCRS)

    # clip raster to buffer sides
    rasterSide1 = raster.rio.clip([bufferSide1])
    rasterSide2 = raster.rio.clip([bufferSide2])
    B = shapely.box(*buffer.bounds)
    rasterTotal = raster.rio.clip([B])

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
        if threshold < rasterSide1.max():
            tree1, coords1        = find_raster_points_threshold(raster1, threshold)
            d1,ind1,  p1, tp, cfs = find_cond_distance(line, confStepSize, tree1, coords1, distanceFactor)
            TP.extend(tp)

        if threshold < rasterSide2.max():
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
    pointsSide1, pointsSide2, TP = [],[], []
    confDistance, confDistance1, confDistance2 = 0,0,0
    lvals = []
    
    for i, l in enumerate(reach.bendLines):

        line_values = extract_slope_along_raster_line(rasterTotal, l, False, '_','_','_',400)
        lineHeight = np.median(line_values[3])
        lvals.append(line_values)


        threshold = ((1/28)*reach.bendWidths[i]*heightFactor) + lineHeight
        D, D1, D2, p1, p2,tp = conf_distance(rasterSide1, rasterSide2, reach.bendWidths[i], threshold, lineHeight, l, stepsize, distanceFactor)
        TP.extend(tp)
        confDistance += D
        confDistance1 += D1
        confDistance2 += D2

        pointsSide1.append(p1)
        pointsSide2.append(p2)

        # plt.plot(np.arange(0, len(line_values[3])),line_values[3])
        # plt.show()
    
    ###################################
    # Plot
    ###################################
    if plot == True:
        f, ax = plt.subplots()
        rasterTotal.plot(ax = ax)
        ax.plot(*line.xy, color = 'black')
        for points in pointsSide1:
            
            for p in points:
                ax.scatter(*p, color = 'red', s= 10)
        for points in pointsSide2:
            for p in points:
                ax.scatter(*p, color = 'darkred', s = 10)
        
        for tp in TP:
            ax.scatter(*tp.xy, s = 5, c = 'green')
        plt.show()
        print(confDistance, line.length)
    

    # adjust to confinement per bend!!!!!!!
    return (confDistance / line.length), (confDistance1 / line.length), (confDistance2 / line.length)