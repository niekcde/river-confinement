"""Compact, full-data map rendering. Rasterize line strokes, never sample bends."""
from pathlib import Path
import base64
from io import BytesIO
import json
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
import plotly.graph_objects as go


EXTENT = (-18050000, 18050000, -9030000, 9030000)
OCEAN_COLOR = '#e4f1f7'
LAND_COLOR = '#d0d0d0'
COAST_COLOR = '#b5b9bb'

def plotly_colorscale(colormap):
    """Use the same Matplotlib palette for Plotly colourbars."""
    from matplotlib.colors import to_hex
    return [[float(t), to_hex(colormap(float(t)))] for t in np.linspace(0, 1, colormap.N)]

def notebook_colormap(path):
    """Standalone scripts read the single COLORMAP definition in the notebook."""
    import ast
    from matplotlib.colors import LinearSegmentedColormap
    notebook = json.loads(Path(path).read_text())
    for cell in notebook['cells']:
        if cell['cell_type'] != 'code':
            continue
        for node in ast.parse(''.join(cell['source'])).body:
            if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'COLORMAP' for t in node.targets):
                return eval(compile(ast.Expression(node.value), str(path), 'eval'),
                            {'LinearSegmentedColormap': LinearSegmentedColormap, 'np': np})
    raise ValueError('Define COLORMAP in the notebook Configuration cell.')

def land_background(land, color=LAND_COLOR, coast_color=COAST_COLOR, name=None):
    """Group Natural Earth polygon exteriors into one separated fill trace."""
    import shapely
    rings = shapely.get_exterior_ring(shapely.get_parts(land.geometry.array))
    coords, ids = shapely.get_coordinates(rings, return_index=True)
    xy = np.full((len(coords) + len(rings), 2), np.nan, dtype=np.float32)
    xy[np.arange(len(coords)) + ids] = coords
    return go.Scatter(x=xy[:,0], y=xy[:,1], mode='lines', fill='toself',
        fillcolor=color, line=dict(color=coast_color, width=.4),
        name=name, hoverinfo='skip', showlegend=False)

def mpl_color(value):
    if value.startswith('rgb('):
        return tuple(float(v)/255 for v in value[4:-1].split(','))
    return value

def hillshade_rgba(path, crs, extent, land_geometries, width=4320, opacity=.35):
    """Warp grayscale relief to the display grid and shade only land.

    Extent uses Matplotlib order (left, right, bottom, top). White relief is
    transparent; darker relief shades the existing land colour. Read through
    a warped VRT so the full-resolution source is never loaded into memory.
    """
    import rasterio
    from rasterio.vrt import WarpedVRT
    from rasterio.enums import Resampling
    from rasterio.transform import from_bounds
    from rasterio.features import rasterize
    if not 0 <= opacity <= 1:
        raise ValueError('Hillshade opacity must be between 0 and 1.')
    left, right, bottom, top = extent
    height = max(1, round(width * (top-bottom)/(right-left)))
    transform = from_bounds(left, bottom, right, top, width, height)
    with rasterio.open(path) as src:
        if src.crs is None or src.count != 1 or src.dtypes[0] != 'uint8':
            raise ValueError('Expected a georeferenced single-band uint8 hillshade.')
        with WarpedVRT(src, crs=crs, transform=transform, width=width,
                       height=height, resampling=Resampling.bilinear,
                       add_alpha=True) as vrt:
            gray = vrt.read(1)
            valid = vrt.dataset_mask() > 0
    land_mask = rasterize(((g, 1) for g in land_geometries if not g.is_empty),
                          out_shape=(height, width), transform=transform,
                          fill=0, dtype='uint8')
    rgba = np.zeros((height, width, 4), dtype=np.uint8)
    rgba[..., 3] = np.rint((255-gray.astype(float))*opacity*valid*land_mask).astype('uint8')
    return rgba


def hillshade_trace(rgba, extent=EXTENT):
    """An image trace sits above polygon fills and below subsequent rivers."""
    from PIL import Image
    buffer = BytesIO()
    Image.fromarray(rgba).save(buffer, format='PNG')
    left, right, bottom, top = extent
    height, width = rgba.shape[:2]
    dx, dy = (right-left)/width, (top-bottom)/height
    return go.Image(source='data:image/png;base64,' + base64.b64encode(buffer.getvalue()).decode(),
                    x0=left+dx/2, y0=top-dy/2, dx=dx, dy=-dy,
                    name='Hillshade', hoverinfo='skip')

def array(value):
    if isinstance(value, dict) and 'bdata' in value:
        return np.frombuffer(base64.b64decode(value['bdata']), dtype=value['dtype'])
    return np.asarray(value, dtype=float)

def read_figure(path):
    text = Path(path).read_text()
    match = list(re.finditer(r'Plotly\.newPlot\(\s*"[^"]+",\s*', text))[-1]
    decoder = json.JSONDecoder()
    data, n = decoder.raw_decode(text[match.end():])
    tail = text[match.end() + n:].lstrip()[1:].lstrip()
    layout, _ = decoder.raw_decode(tail)
    return go.Figure(data=data, layout=layout)

def compact_map(fig, png_path, width=4320):
    """Draw every river segment with round caps; leave land/legends as Plotly vectors.

    Fixed 4320 px world image gives ~8.4 km/pixel at the equator. Zooming the
    compact map magnifies this raster; use *_vectors.html for geometry-level zoom.
    """
    data = fig.to_plotly_json()['data']
    rivers = [t for t in data if t.get('mode') == 'lines' and t.get('fill') != 'toself']
    background = [t for t in data if t.get('fill') == 'toself']
    relief = [t for t in data if t.get('type') == 'image']
    bars = [t for t in data if t.get('marker', {}).get('showscale')]
    x0, x1, y0, y1 = EXTENT
    dpi = 180
    canvas = plt.figure(figsize=(width/dpi, width/2/dpi), dpi=dpi)
    ax = canvas.add_axes([0, 0, 1, 1])
    ax.set(xlim=(x0,x1), ylim=(y0,y1)); ax.set_axis_off()
    segment_count = 0
    for t in rivers:
        x, y = array(t['x']), array(t['y'])
        finite = np.isfinite(x) & np.isfinite(y)
        segment_count += int(np.count_nonzero(finite & np.r_[True, ~finite[:-1]]))
        ax.plot(x, y, color=mpl_color(t['line']['color']), linewidth=.55,
                solid_capstyle='round', antialiased=True)
    buffer = BytesIO()
    canvas.savefig(buffer, format='png', transparent=True, dpi=dpi)
    plt.close(canvas)
    raster_bytes = buffer.getvalue()
    uri = 'data:image/png;base64,' + base64.b64encode(raster_bytes).decode()
    compact = go.Figure(layout=fig.layout)
    for t in background:
        compact.add_trace(go.Scatter(**{k:v for k,v in t.items() if k != 'type'}))
    for t in relief:
        compact.add_trace(go.Image(**{k:v for k,v in t.items() if k != 'type'}))
    # Raster contains every bend, in exactly the same layer order as the vector map.
    compact.add_layout_image(dict(source=uri, xref='x', yref='y', x=x0, y=y1,
        sizex=x1-x0, sizey=y1-y0, sizing='stretch', layer='above'))
    for t in rivers:
        if t.get('showlegend', True):
            compact.add_trace(go.Scatter(x=[None], y=[None], mode='lines', name=t['name'],
                line=dict(color=t['line']['color'], width=2), showlegend=True))
    for t in bars:
        compact.add_trace(go.Scatter(**{k:v for k,v in t.items() if k != 'type'}))
    compact.update_layout(meta=dict(rendering='All bends rasterized; no sampling',
        river_line_parts=segment_count, raster_width_px=width))
    compact.update_layout(legend=dict(itemclick=False, itemdoubleclick=False))
    # Browser-independent static export, including ALL river strokes.
    static, ax = plt.subplots(figsize=(16, 8.8), dpi=220)
    page_color = fig.layout.paper_bgcolor or 'white'
    static.set_facecolor(page_color)
    static.subplots_adjust(left=.01, right=.9 if bars else .99, bottom=.055, top=.91)
    ax.set(xlim=(x0,x1), ylim=(y0,y1), aspect='equal'); ax.set_axis_off()
    for t in background:
        ax.fill(array(t['x']), array(t['y']), facecolor=mpl_color(t['fillcolor']),
                edgecolor=mpl_color(t.get('line', {}).get('color', t['fillcolor'])), linewidth=.25)
    raster = plt.imread(BytesIO(raster_bytes), format='png')
    for t in relief:
        relief_image = plt.imread(BytesIO(base64.b64decode(t['source'].split(',', 1)[1])), format='png')
        h, w = relief_image.shape[:2]
        left, top = t['x0']-t['dx']/2, t['y0']-t['dy']/2
        ax.imshow(relief_image, extent=(left, left+w*t['dx'], top+h*t['dy'], top),
                  origin='upper', interpolation='nearest', zorder=2)
    ax.imshow(raster, extent=EXTENT, origin='upper', interpolation='nearest', zorder=3)
    labels = [t for t in rivers if t.get('showlegend', True)]
    if fig.layout.legend.traceorder == 'reversed': labels.reverse()
    if labels:
        ax.legend(handles=[Line2D([0],[0], color=mpl_color(t['line']['color']), lw=2, label=t['name']) for t in labels],
            loc='lower left', frameon=True, fontsize=10)
    if bars:
        marker = bars[0]['marker']; spec = marker['colorbar']
        cax = static.add_axes([.925,.19,.016,.61])
        from matplotlib.colors import LinearSegmentedColormap
        cmap = LinearSegmentedColormap.from_list('map_colors',
            [(position, mpl_color(color)) for position, color in marker['colorscale']])
        cb = static.colorbar(plt.cm.ScalarMappable(norm=Normalize(marker['cmin'],marker['cmax']),
                                                   cmap=cmap),cax=cax)
        cb.set_ticks(spec['tickvals'], labels=spec['ticktext'])
        cax.set_title('Mean ratio',fontsize=10)
    static.suptitle(fig.layout.title.text, fontsize=16, x=.04, ha='left')
    # static.text(.02,.018,'Mollweide • Natural Earth 1:10m land • ALL bends • mean = (inner + outer) / 2', fontsize=9)
    static.savefig(png_path, dpi=500, facecolor=page_color)
    plt.close(static)
    print(f'{Path(png_path).name}: {segment_count:,} line parts, no sampling', flush=True)
    return compact

if __name__ == '__main__':
    import geopandas as gpd
    root=Path(__file__).resolve().parents[1]/'results/global_confinement'
    land=gpd.read_file(f"zip://{Path.home() / 'Downloads/10m_physical.zip'}!ne_10m_land.shp").to_crs(
        '+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs')
    background=land_background(land)
    for name in ['mollweide_classes','mollweide_log','mollweide_log1p']:
        original=root/f'{name}_vectors.html'
        current=root/f'{name}.html'
        if not original.exists(): current.rename(original)
        fig=read_figure(original)
        fig=go.Figure(data=[background]+[t for t in fig.data if t.fill != 'toself'], layout=fig.layout)
        fig.update_layout(paper_bgcolor=OCEAN_COLOR, plot_bgcolor=OCEAN_COLOR)
        for annotation in fig.layout.annotations:
            annotation.text=annotation.text.replace('1:110m', '1:10m')
        fig.write_html(original,include_plotlyjs=True)
        result=compact_map(fig,root/f'{name}.png')
        result.write_html(current,include_plotlyjs=True)
        result.write_json(root/f'{name}_compact.json')
        print('Saved',current,flush=True)
