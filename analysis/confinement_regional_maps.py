"""Equal-scale UTM region zooms from ORIGINAL confinement bend geometries.

Optional validated bend-ID lookups restrict chosen windows to a catchment;
unfiltered windows can include neighbouring basins. No bend sampling or geometry
simplification is applied. The notebook supplies one colormap for every map.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
from pyproj import CRS, Transformer, Proj
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
from PIL import Image

REGIONS = [
    dict(name='Rhine', slug='rhine', lon=8.0, lat=50.0, epsg=32632),
    dict(name='Amazon', slug='amazon', lon=-60.5, lat=-3.5, epsg=32720),
    dict(name='Ob', slug='ob', lon=73.5, lat=61.0, epsg=32643),
    dict(name='Brahmaputra', slug='brahmaputra', lon=93.0, lat=27.0, epsg=32646),
]
LINEWIDTH = .42
def rgb(color):
    return tuple(float(v)/255 for v in color[4:-1].split(',')) if color.startswith('rgb(') else color


def polygons(ax, geometries, color, edge, zorder):
    """Draw both exterior and interior rings, preserving polygon holes."""
    for polygon in shapely.get_parts(geometries):
        if polygon.is_empty or polygon.geom_type != 'Polygon':
            continue
        polygon = shapely.geometry.polygon.orient(polygon, sign=1.0)
        paths = []
        for ring in [polygon.exterior, *polygon.interiors]:
            xy = np.asarray(ring.coords)
            codes = np.full(len(xy), MplPath.LINETO, dtype=np.uint8)
            codes[0] = MplPath.MOVETO; codes[-1] = MplPath.CLOSEPOLY
            paths.append(MplPath(xy, codes))
        patch = PathPatch(MplPath.make_compound_path(*paths), facecolor=color,
                          edgecolor=edge, linewidth=.25, zorder=zorder)
        ax.add_patch(patch)


def clipped_background(path, bbox, crs, window):
    frame = gpd.read_file(path, bbox=bbox)
    # Clip geographic background BEFORE projection: distant polygons can fail
    # in a local UTM zone. Add a small geographic buffer before exact UTM clip.
    geo_box = shapely.box(bbox[0]-.2, bbox[1]-.2, bbox[2]+.2, bbox[3]+.2)
    frame.geometry = shapely.intersection(shapely.make_valid(frame.geometry.array), geo_box)
    # Geographic clipping creates long straight edges. Densify before projecting
    # so latitude/longitude edges follow their UTM curves, rather than chords
    # that cut across the map (especially in large, high-latitude windows).
    frame.geometry = shapely.segmentize(frame.geometry.array, max_segment_length=.1)
    frame = frame[~frame.geometry.is_empty].to_crs(crs)
    frame.geometry = shapely.intersection(shapely.make_valid(frame.geometry.array), window)
    return frame.geometry.array


def render_regions(source, land_dir, out, regions=REGIONS, REGIONAL_WINDOW_KM=1000,
                   log_min=.5, cap=50., ocean_color='#9ebcd8',
                   land_color='#dfd3d1', coast_color='#dfd3d1',
                   page_color='white', nan_color=None, dpi=500, *, colormap,
                   catchment_lookup=None, catchments=None, render_only=None,
                   hillshade_path=None, hillshade_opacity=.35, linewidths=LINEWIDTH):
    """Render maps with river widths in points.

    linewidths accepts one number, a list in regions order, or a dictionary
    keyed by region slug (useful when rendering only one selected region).
    """
    if isinstance(linewidths, dict):
        missing = [r['slug'] for r in regions if r['slug'] not in linewidths]
        if missing:
            raise ValueError(f'Missing linewidths for regions: {missing}')
        region_widths = np.asarray([linewidths[r['slug']] for r in regions], dtype=float)
    elif np.isscalar(linewidths):
        region_widths = np.full(len(regions), float(linewidths))
    else:
        region_widths = np.asarray(linewidths, dtype=float)
    if region_widths.shape != (len(regions),):
        raise ValueError('Provide one linewidth per region, in regions order.')
    if not np.all(np.isfinite(region_widths) & (region_widths > 0)):
        raise ValueError('Linewidths must be finite positive numbers (points).')
    source, land_dir, out = Path(source), Path(land_dir), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    assert 0 < log_min < cap
    assert isinstance(REGIONAL_WINDOW_KM, (int, float)) or isinstance(REGIONAL_WINDOW_KM, list)
    limits = np.log10([log_min, cap])
    cmap = colormap
    if nan_color is None:
        nan_color = cmap(1.0)
    catchments = catchments or {}
    lookup = None
    if catchments:
        if catchment_lookup is None:
            raise ValueError('Provide a validated bend catchment lookup for catchment filtering.')
        lookup = pd.read_csv(catchment_lookup,
            usecols=['file','bendID','combined_reach_id','catchment','assignment_status'],
            dtype={'file':'string','bendID':'string','combined_reach_id':'int64',
                   'catchment':'string','assignment_status':'string'})
        assert not lookup.duplicated(['file','bendID']).any()
        lookup = lookup.rename(columns={'combined_reach_id':'lookup_combined_reach_id'})
    previous_path = out / 'regional_windows.json'
    previous = {r['slug']:r for r in json.loads(previous_path.read_text())} if previous_path.exists() else {}
    metadata=[]; image_paths=[]
    c = 0
    for region in regions:
        river_linewidth = float(region_widths[c])
        window_km = REGIONAL_WINDOW_KM[c] if isinstance(REGIONAL_WINDOW_KM, list) else REGIONAL_WINDOW_KM
        c += 1
        if render_only is not None and region['slug'] not in render_only:
            path = out / f'{region["slug"]}_log_utm.png'
            if not path.exists() or region['slug'] not in previous:
                raise ValueError(f'No existing output to preserve for {region["slug"]}.')
            metadata.append(previous[region['slug']]); image_paths.append(path)
            print(f'Preserved existing {region["name"]} image',flush=True)
            continue
        crs=CRS(region['epsg'])
        forward=Transformer.from_crs(4326,crs,always_xy=True)
        inverse=Transformer.from_crs(crs,4326,always_xy=True)
        cx,cy=forward.transform(region['lon'],region['lat'])
        half=window_km*1000/2
        extent=(cx-half,cy-half,cx+half,cy+half)
        bbox=inverse.transform_bounds(*extent,densify_pts=81)
        window=shapely.box(*extent)
        print(f"{region['name']}: reading original bends in {crs.name}",flush=True)
        target = catchments.get(region['slug'])
        columns=['bendID','ER_inn','ER_out']
        if target:
            columns += ['file','combined_reach_id_x']
        bends=gpd.read_file(source,layer='global_clusters_bend',bbox=bbox,columns=columns)
        bends=bends.to_crs(crs)
        bends=bends.loc[bends.geometry.intersects(window)].copy()
        bends.geometry=shapely.intersection(bends.geometry.array,window)
        bends=bends.loc[(~bends.geometry.is_empty) & (bends.geometry.length>0)].copy()
        unfiltered_count = len(bends)
        unassigned_count = 0
        if target:
            bends = bends.merge(lookup,on=['file','bendID'],how='left',validate='many_to_one')
            if bends.assignment_status.isna().any():
                raise ValueError('Some bend IDs are absent from the lookup; rebuild it for this GeoPackage.')
            assert bends.combined_reach_id_x.eq(bends.lookup_combined_reach_id).all()
            unassigned_count = int(bends.assignment_status.ne('unambiguous').sum())
            bends = bends.loc[bends.catchment.eq(target)].copy()
            assert bends.assignment_status.eq('unambiguous').all()
            bends[['file','bendID','lookup_combined_reach_id','catchment']].rename(
                columns={'lookup_combined_reach_id':'combined_reach_id'}).to_csv(
                    out/f'{region["slug"]}_selected_bend_ids.csv.gz',index=False)
        assert len(bends)>0
        ratio=((bends.ER_inn+bends.ER_out)/2).to_numpy()
        finite=np.isfinite(ratio)
        assert not np.any(finite & (ratio<0))
        assert not np.any(~finite & ~np.isnan(ratio))
        positive=finite & (ratio>0)
        normalized=np.zeros(len(ratio))
        normalized[positive]=(np.log10(np.clip(ratio[positive],log_min,cap))-limits[0])/(limits[1]-limits[0])
        bands=np.floor(normalized*96).clip(0,95).astype(int)
        # Split only multipart lines created by clipping, never join bends.
        parts, parents=shapely.get_parts(bends.geometry.array,return_index=True)
        line_mask=shapely.get_type_id(parts)==1
        parts,parents=parts[line_mask],parents[line_mask]
        assert len(np.unique(parents))==len(bends), 'A visible bend lost its line geometry.'
        segments=[np.asarray(g.coords) for g in parts]
        land=clipped_background(land_dir/'ne_10m_land.shp',bbox,crs,window)
        ocean=clipped_background(land_dir/'ne_10m_ocean.shp',bbox,crs,window)
        fig=plt.figure(figsize=(8,8.8),dpi=dpi,facecolor=page_color)
        # Exactly 6.24 x 6.24 inches for every 1000 x 1000 km map frame.
        ax=fig.add_axes([.045,.15,.78,.78*8/8.8])
        ax.set(xlim=(extent[0],extent[2]),ylim=(extent[1],extent[3]),aspect='equal')
        ax.set_facecolor(page_color);ax.set_xticks([]);ax.set_yticks([])
        for spine in ax.spines.values():spine.set_color('#9c9c9c');spine.set_linewidth(.4)
        polygons(ax,ocean,ocean_color,ocean_color,0)
        polygons(ax,land,land_color,coast_color,1)
        if hillshade_path is not None:
            from confinement_rendering import hillshade_rgba
            image_extent = (extent[0], extent[2], extent[1], extent[3])
            relief = hillshade_rgba(hillshade_path, crs, image_extent, land,
                                   width=round(6.24*dpi), opacity=hillshade_opacity)
            ax.imshow(relief, extent=image_extent, origin='upper', zorder=1.5)
        painted=0
        def draw(mask,color):
            nonlocal painted
            selected=np.flatnonzero(mask[parents])
            if len(selected):
                collection=LineCollection([segments[i] for i in selected],colors=color,
                    linewidths=river_linewidth,capstyle='round',zorder=2)
                ax.add_collection(collection);painted+=len(selected)
        draw(np.isnan(ratio),nan_color)
        draw(finite & (ratio==0),'#888888')
        for band in reversed(range(96)):
            draw(positive & (bands==band),cmap((band+.5)/96))
        assert painted==len(parts)
        bar_km=100 if window_km>=400 else window_km/5
        xb,yb=extent[0]+.06*2*half,extent[1]+.055*2*half
        ax.plot([xb,xb+bar_km*1000],[yb,yb],color='black',lw=2,zorder=4)
        ax.text(xb+bar_km*500,yb+.018*2*half,f'{bar_km:g} km',ha='center',fontsize=9,
                bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=1),zorder=4)
        ax.annotate('Grid N',xy=(.95,.95),xytext=(.95,.86),xycoords='axes fraction',
            ha='center',fontsize=8,arrowprops=dict(arrowstyle='-|>',color='#222'),zorder=4)
        cax=fig.add_axes([.865,.24,.024,.56])
        cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(*limits),cmap=cmap),cax=cax)
        ticks=[v for v in [.5,1,2,5,10,15,25,50] if log_min<=v<=cap]
        cb.set_ticks(np.log10(ticks),labels=[f'{v:g}' for v in ticks])
        cax.set_title('Mean ratio',fontsize=9,pad=10)
        cb.ax.tick_params(labelsize=9)
        fig.text(.045,.955,region['name'],fontsize=17,weight='semibold')
        fig.text(.045,.922,f'{crs.name} · EPSG:{region["epsg"]}',fontsize=10)
        # fig.text(.045,.103,f'{window_km:g} × {window_km:g} km · log colours {log_min:g}–{cap:g} · NaN = upper-end colour',fontsize=9)
        selection_note = ('HydroBASINS catchment only; cropped to the common window.' if target
                          else 'Regional window; neighbouring catchments may appear.')
        fig.text(.045,.070,selection_note,fontsize=9,color='#555')
        png=out/f'{region["slug"]}_log_utm.png'
        fig.savefig(png,dpi=dpi,facecolor=page_color)
        plt.close(fig);image_paths.append(png)
        # Quantify local distortion when windows extend beyond a zone boundary.
        gx,gy=np.meshgrid(np.linspace(extent[0],extent[2],11),np.linspace(extent[1],extent[3],11))
        lons,lats=inverse.transform(gx.ravel(),gy.ravel())
        factors=np.asarray(Proj(crs).get_factors(lons,lats).meridional_scale)
        row={**region,'window_km':window_km,'log_min':log_min,'log_max':cap,
             'river_linewidth_pt':river_linewidth,
             'colormap':cmap.name,
             'catchment_filter':target,'unfiltered_bends':unfiltered_count,
             'excluded_unassigned_bends':unassigned_count,
             'visible_bends':len(bends),'drawn_line_parts':len(parts),
             'missing_mean':int(np.isnan(ratio).sum()),'zero_mean':int((ratio==0).sum()),
             'x_min':extent[0],'x_max':extent[2],'y_min':extent[1],'y_max':extent[3],
             'point_scale_min':float(factors.min()),'point_scale_max':float(factors.max()),
             'source_vertices_clipped':int(shapely.get_num_coordinates(parts).sum()),
             'png':str(png),'frame_width_inches':6.24,'dpi':dpi}
        metadata.append(row)
        print(f"Saved {png.name}: {len(bends):,} bends; point scale {factors.min():.6f}–{factors.max():.6f}",flush=True)
    # A compact contact sheet with equal pixel dimensions for each panel.
    thumbs=[]
    for path in image_paths:
        with Image.open(path) as im:
            thumbs.append(im.convert('RGB').resize((1400,1540),Image.Resampling.LANCZOS))
    sheet=Image.new('RGB',(2800,3080),'white')
    for i,im in enumerate(thumbs):sheet.paste(im,((i%2)*1400,(i//2)*1540))
    sheet.save(out/'regional_log_utm_comparison.png')
    pd.DataFrame(metadata).to_csv(out/'regional_windows.csv',index=False)
    (out/'regional_windows.json').write_text(json.dumps(metadata,indent=2))
    return pd.DataFrame(metadata)
