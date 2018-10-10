import numpy as np
import os
import sys

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import prepare_EOlab_layers as EO

sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/PotentialBiomass/src')
import geospatial_utility_tools as geo
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/')
import auxilliary_functions as aux

# Get perceptionally uniform colourmaps
sys.path.append('/home/dmilodow/DataStore_DTM/FOREST2020/EOdata/EO_data_processing/src/plot_EO_data/colormap/')
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='inferno', cmap=cmaps.inferno)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
plt.register_cmap(name='magma', cmap=cmaps.magma)
plt.set_cmap(cmaps.viridis)


import shapefile
sf = shapefile.Reader("/home/dmilodow/DataStore_DTM/EOlaboratory/Areas/NaturalEarth/10m_cultural/ne_10m_admin_0_countries")

vertices = []
for shape_rec in sf.shapeRecords():
    if shape_rec.record[3] == 'Ghana':
        pts = shape_rec.shape.points
        prt = list(shape_rec.shape.parts) + [len(pts)]
        for i in range(len(prt) - 1):
            vertices.append([])
            for j in range(prt[i], prt[i+1]):
                vertices[i].append((pts[j][0], pts[j][1]))


plt.figure(1, facecolor='White',figsize=[2, 1])
plt.show()

DATADIR = '/home/dmilodow/DataStore_DTM/FOREST2020/PotentialBiomassV2/output/'
SAVEDIR = '/home/dmilodow/DataStore_DTM/EOlaboratory/EOlab/GhanaPotentialAGB/'
NetCDF_file = 'EAFR_v002_AGBpot_mean_WC2_SOILGRIDS_GridSearch.nc'

ds,geoTrans = EO.load_NetCDF(DATADIR+NetCDF_file,lat_var = 'lat', lon_var = 'lon')
lat = np.asarray(ds.variables['lat'])
lon = np.asarray(ds.variables['lon'])
vars = ['AGB_mean','AGBpot_mean','training', 'areas']

# create clipping mask
S=4.5
N=11.5
E= 1.5
W=-3.5
lat_mask = np.all((lat<=N,lat>=S),axis=0)
lon_mask = np.all((lon<=E,lon>=W),axis=0)
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

# now loop through all polygons, and find the pixels which fall within them
country_mask = np.zeros((n_lat,n_lon))

for vv in range(0,len(vertices)):
      temp1,temp2,inside = aux.points_in_poly(longrid.ravel(),latgrid.ravel(),vertices[vv])
      country_mask[inside.reshape(n_lat,n_lon)] = 1
"""
for ilat in range(0,lat_clip.size):
  print ilat, lat_clip.size
  for ilon in range(0,lon_clip.size):
    shape = 0
    while shape<len(vertices):

      if aux.point_in_poly(lon_clip[ilon],lat_clip[ilat],vertices[shape]):
        country_mask[ilat,ilon]=1
        shape = len(vertices)
      shape+=1
"""

# apply clip
dataset = {}
for vv in range(0,len(vars)):
    dataset[vars[vv]]= np.asarray(ds.variables[vars[vv]])[array_mask]
    dataset[vars[vv]][country_mask==0]=-9999

# create new geoTrans object to account for clipping
geoTrans[0] = np.min(lon_clip)
if geoTrans[5]>0:
  geoTrans[3]=np.min(lat_clip)
else:
  geoTrans[3]=np.max(lat_clip)

# sequestration potential is defined by pixels with positive potential biomass that
# are not already forests
dataset['seqpot_mean'] = dataset['AGBpot_mean']-dataset['AGB_mean']
#dataset['seqpot_mean'][dataset['training']==1] = 0.
dataset['seqpot_mean'][dataset['seqpot_mean']<0] = 0.
dataset['seqpot_mean'][dataset['AGB_mean']==-9999] = -9999.

dataset['training'][dataset['training']<1] = -9999.

vars = ['AGB_mean','AGBpot_mean','seqpot_mean','training']
cmaps = ['viridis','viridis','plasma','viridis']
ulims = [400.,400.,200.,1.]
llims = [0.,0.,0.,0.]
axis_labels = ['AGB$_{obs}$ / Mg(C) ha$^{-1}$', 'AGB$_{potential}$ / Mg(C) ha$^{-1}$', 'Sequestration potential / Mg(C) ha$^{-1}$', 'Training sample']

for vv in range(0,len(vars)):
    print( vars[vv])
    file_prefix = SAVEDIR + 'ghana_' + vars[vv]

    # delete existing dataset if present
    if 'ghana_'+vars[vv]+'_data.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % (SAVEDIR+'ghana_'+vars[vv]+'_data.tif'))
    if 'ghana_'+vars[vv]+'_display.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % (SAVEDIR+'ghana_'+vars[vv]+'_display.tif'))

    EO.write_array_to_display_layer_GeoTiff(dataset[vars[vv]], geoTrans, file_prefix, cmaps[vv], ulims[vv], llims[vv])
    if vars!='training':
        EO.plot_legend(cmaps[vv],ulims[vv],llims[vv],axis_labels[vv], file_prefix,extend='max')

"""
rows, cols = dataset[vars[0]].shape
latitude = np.arange(geoTrans_rs[3],rows*geoTrans_rs[5]+geoTrans_rs[3],geoTrans_rs[5])[:rows]
longitude =  np.arange(geoTrans_rs[0],cols*geoTrans_rs[1]+geoTrans_rs[0],geoTrans_rs[1])[:cols]
areas = geo.calculate_cell_area_array(latitude,longitude, area_scalar = 1./10.**4,cell_centred=False)
"""

# loop through the variables, multiplying by cell areas to give values in Mg
for vv in range(0,len(vars)):
    print(vars[vv])

    if 'ghana_'+vars[vv]+'_total_data.tif' in os.listdir(SAVEDIR):
        os.system("rm %s" % (SAVEDIR+'ghana_'+vars[vv]+'_total_data.tif'))

    file_prefix = SAVEDIR + 'ghana_' + vars[vv] + '_total'

    out_array = dataset[vars[vv]] * dataset['areas']
    out_array[dataset[vars[vv]]==-9999]=-9999
    EO.write_array_to_data_layer_GeoTiff(out_array, geoTrans_rs, file_prefix)
    out_array=None

# Also want to write cell areas to file.  However, as this will be compared against other layers, need to carry across
# nodata values
areas_out = dataset['areas']
areas_out[np.asarray(dataset[vars[0]])==-9999] = -9999
if 'ghana_cell_areas_data.tif' in os.listdir(SAVEDIR):
    os.system("rm %s" % (SAVEDIR+'ghana_cell_areas_data.tif'))
area_file_prefix = SAVEDIR + 'ghana_cell_areas'
EO.write_array_to_data_layer_GeoTiff(areas_out, geoTrans_rs, area_file_prefix)
