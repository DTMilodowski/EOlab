"""
PotentialAGB_app.py
This function creates the required layers for potential AGB EO lab applications.
Specific layers produced:
    - Observed AGB
    - Potential AGB
    - Deficit
    - Sequestration Potential
    - Grid cell areas (for on-the-fly observations)
Both data (WGS84) and display (RGB, Web Mercator) layers are produced, with an
associated legend for the latter.

Input layers are assumed to be netcdf files produced with the scripts in the
PotentialBiomassRFR library.

Input arguments help to point to the correct files for the country of interest.
If the "country" option is specified as "None", then by default no clipping of
the data will be undertaken, otherwise the data will be clipped to the national
boundaries of the specified country.
"""
## READ INPUT ARGUMENTS
from  input_arguments_mexico import country,outfileID, DATADIR, SAVEDIR, NetCDF_file, ne_shp,\
                                    Cscale, map_vars,cmaps,ulims,llims,axis_labels

## IMPORT PACKAGES
import numpy as np
import pandas as pd
import os
import sys

# Plotting libraries
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import seaborn as sns
plt.cla()

# GIS libraries
import xarray as xr
import geopandas
from geopandas.tools import sjoin
from cartopy.io.shapereader import Reader
from shapely.geometry import Point

# custom libraries
sys.path.append('../')
import prepare_EOlab_layers as EO

# create colormap for training data
colours = np.asarray(['#46E900','#1A2BCE'])
training_labels = np.asarray(['initial','supplementary'])
training_cmap = ListedColormap(sns.color_palette(colours).as_hex())
plt.register_cmap(name='training_cmap', cmap=training_cmap)

## READ IN DATA
ds = xr.open_dataset('%s%s' % (DATADIR,NetCDF_file))
lat = ds.variables['lat'].values
lon = ds.variables['lon'].values
dlat = lat[1]-lat[0]
dlon = lon[1]-lon[0]

vars = ['AGBobs']
vars.append(list(ds.data_vars)[-2])
vars.append(list(ds.data_vars)[-1])

## CLIP IF REQUIRED
if country!="None":
    # get bounds of national border shapefile to clip raster to required extent
    sf = Reader(ne_shp)
    boundaries = [i for i in sf.records() if i.attributes['admin'] == country]
    N=0.;S=0.;W=0.;E=0
    count=0
    for b in boundaries:
        if count == 0:
            N = b.bounds[3];S = b.bounds[1];E = b.bounds[2];W = b.bounds[0]
            count+=1
        else:
            N=max(N,b.bounds[3]);E=max(E,b.bounds[2]);W=min(W,b.bounds[0]);S=min(S,b.bounds[1])

    # clip raster based on this extent
    lat_mask = np.all((lat<=N+abs(dlat),lat>=S-abs(dlat)),axis=0)
    lon_mask = np.all((lon<=E+abs(dlon),lon>=W-abs(dlon)),axis=0)
    n_lat = lat_mask.sum()
    n_lon = lon_mask.sum()
    lat_clip = lat[lat_mask]
    lon_clip = lon[lon_mask]
    array_mask = np.ix_(lat_mask,lon_mask)

    latgrid = np.zeros((n_lat,n_lon))
    longrid = np.zeros((n_lat,n_lon))
    for ilat in range(0,n_lat):
      longrid[ilat,:] = lon_clip.copy()

    for ilon in range(0,n_lon):
      latgrid[:,ilon] = lat_clip.copy()

    # now mask pixels that fall outside of country boundary
    inside = np.zeros(n_lat*n_lon)
    testlon = longrid.reshape(n_lat*n_lon)
    testlat = latgrid.reshape(n_lat*n_lon)

    boundaries = geopandas.GeoDataFrame.from_file(ne_shp)
    boundaries = boundaries[boundaries['admin'] == country]

    df = pd.DataFrame({'latitude':testlat,'longitude':testlon})
    df['Coordinates'] = list(zip(df.longitude, df.latitude))
    df['Coordinates'] = df['Coordinates'].apply(Point)
    point = geopandas.GeoDataFrame(df, geometry='Coordinates')
    pointInPolys = sjoin(point, boundaries, how='left')
    country_mask = np.reshape(np.asarray(pointInPolys.admin)==country,(n_lat,n_lon))

    # Save into dictionary of processed layers
    dataset = {}
    for vv in range(0,len(vars)):
        dataset[vars[vv]]= np.asarray(ds.variables[vars[vv]])[array_mask]
        dataset[vars[vv]][country_mask==0]=np.nan
        dataset[vars[vv]][dataset[vars[vv]]==-9999]=np.nan

    dataset['AGBpot'] = dataset.pop(vars[1])
    dataset['training'] = dataset.pop(vars[2])

    # create new geoTrans object to account for clipping
    geoTrans = [0,lon_clip[1]-lon_clip[0],0,0,0,lat_clip[1]-lat_clip[0]]
    geoTrans[0] = np.min(lon_clip)
    if geoTrans[5]>0:
      geoTrans[3]=np.min(lat_clip)
    else:
      geoTrans[3]=np.max(lat_clip)

## CALCULATE ADDITIONAL LAYERS
# sequestration potential is defined by pixels with positive potential biomass that
# are not already forests
dataset['AGBobs']*=Cscale
dataset['AGBpot']*=Cscale
dataset['AGBdef'] = dataset['AGBpot']-dataset['AGBobs']
dataset['AGBseq'] = dataset['AGBdef'].copy()
dataset['AGBseq'][dataset['AGBseq']<0] = 0.
dataset['training'][dataset['training']<1] = np.nan

# calculate cell areas
areas = np.zeros((n_lat,n_lon))
res = np.abs(lat_clip[1]-lat_clip[0])
for la,latval in enumerate(lat_clip):
    areas[la]= (6371e3)**2 * ( np.deg2rad(0+res/2.)-np.deg2rad(0-res/2.) ) * (np.sin(np.deg2rad(latval+res/2.))-np.sin(np.deg2rad(latval-res/2.)))
areas[np.isnan(dataset['AGBobs'])] = np.nan
areas/=10.**4 # m2 to ha

for vv in range(0,len(map_vars)):
    print(map_vars[vv])
    file_prefix = SAVEDIR + outfileID + '_' + map_vars[vv]

    # remove existing datsets
    if (file_prefix+'data.tif') in os.listdir(SAVEDIR):
        os.system("rm %s" % (file_prefix+'_data.tif'))
    if (file_prefix+'data.tif') in os.listdir(SAVEDIR):
        os.system("rm %s" % (file_prefix+'_display.tif'))

    out_array = dataset[map_vars[vv]].copy()
    out_array[np.isnan(out_array)]=-9999.
    EO.write_array_to_display_layer_GeoTiff(out_array, geoTrans, file_prefix, cmaps[vv], ulims[vv], llims[vv])
    plt.cla()
    if map_vars[vv]=='training':
        EO.plot_legend_listed(training_cmap,training_labels,axis_labels[vv], file_prefix)
    elif map_vars[vv]=='AGBdef':
        EO.plot_legend(cmaps[vv],ulims[vv],llims[vv],axis_labels[vv], file_prefix,extend='both')
    else:
        EO.plot_legend(cmaps[vv],ulims[vv],llims[vv],axis_labels[vv], file_prefix,extend='max')

    # multiply by cell areas to give values in Mg
    if map_vars[vv]!='training':
        file_prefix = SAVEDIR + outfileID + '_' + map_vars[vv] + '_total'
        if (file_prefix+'_data.tif') in os.listdir(SAVEDIR):
            os.system("rm %s" % (file_prefix+'_data.tif'))

        out_array = dataset[map_vars[vv]] * areas
        out_array[np.isnan(out_array)]=-9999
        EO.write_array_to_data_layer_GeoTiff(out_array, geoTrans, file_prefix)
        out_array=None

# Also want to write cell areas to file.  However, as this will be compared against other layers, need to carry across
# nodata values
if (outfileID+'_cell_areas_data.tif') in os.listdir(SAVEDIR):
    os.system("rm %s" % (SAVEDIR+outfileID + '_cell_areas_data.tif'))
area_file_prefix = SAVEDIR + outfileID + '_cell_areas'
EO.write_array_to_data_layer_GeoTiff(areas, geoTrans, area_file_prefix)
